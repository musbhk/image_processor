classdef imageprocessor
    % IMAGEPROCESSOR 
    % upload image of starfield and extracts centroids
    %
    % typical use:
    %   ip = imageprocessor();
    %   img = ip.load_image('stellarium_orion.png');
    %   stars = ip.detect_stars(img);   % [N x 3]: u, v, brightness
    %
    % stars(:,1) = horizontal coord (column) in pixel
    % stars(:,2) = vertical coord (row) in pixel
    % stars(:,3) = estimated brightness (somma intensità nel blob)

    properties
        threshold_sigma    % how many std's above avg
        min_blob_size      %min. pixels to consider a blob a star
        max_blob_size      % max. pixels (to exclude planets/Moon)
        saturation_value   % saturation value (0-255 RGB uint8)
    end

    methods
        function obj = imageprocessor()
            % class constructor w/ default values 
            obj.threshold_sigma  = 3.5;
            obj.min_blob_size    = 3;
            obj.max_blob_size    = 50;
            obj.saturation_value = 255;
        end

        function img = load_image(obj, filename)
            % uploads image and converts in grayscale double [0, 1].
            %
            % Step:
            %  1. imread to read PNG/JPG
            %  2. if RGB, converts in grayscale (luminance)
            %  3. converts in double and normalizes to [0, 1]

            raw = imread(filename);

            if size(raw, 3) == 3
                % RGB → perceived luminance (ITU-R BT.601 formula)
                R = double(raw(:,:,1));
                G = double(raw(:,:,2));
                B = double(raw(:,:,3));
                gray = 0.299*R + 0.587*G + 0.114*B;
            else
                gray = double(raw);
            end

            % normalizes in [0, 1]
            img = gray / 255.0;     % this represents a 2D matrix of luminosity, in which every number is the luminosity of a pixel in [0,1]
        end

        function stars = detect_stars(obj, img)
            % detects star centroids w/ threshold + connected components.
            %
            % output: [N x 3] = [u, v, brightness]

            % === step 1: global threshold ===
            mu  = mean(img(:));  % this is the mean luminosity of the whole image. in astronomical images, the dark sky covers most of the pixels, so this value actually shows us the background noise of the image (camera+zodiacal light+unresolved stars)
            sig = std(img(:));   % this is the std of the background noise. if a pixel is off by 1*sig from mu, it could be noise, but if it's 5*sig off from mu, then it certainly is something real
            T   = mu + obj.threshold_sigma * sig; % this combines these two measures of noise to give us a threshold of noise. any pixel with brghtness > T is most certainly a star

            mask = img > T;   % logic matrix: "star" pixel = true

            % === step 2: connected components ===
            % finds blobs (groups of adjacent pixels that are "true")
            % simple implementation: flood fill o bwlabel-like.
            [labels, n_blobs] = obj.label_connected_components(mask);

            % === step 3: centroid of each blob ===
            stars = zeros(n_blobs, 3);
            valid = false(n_blobs, 1);

            for k = 1:n_blobs
                [rows, cols] = find(labels == k);
                blob_size = numel(rows);

                % filters blob dimensions
                if blob_size < obj.min_blob_size || blob_size > obj.max_blob_size
                    continue;
                end

                % intensity of blob's pixels
                idx = sub2ind(size(img), rows, cols);
                intensities = img(idx);

                total = sum(intensities);

                % weighed centroid (sub-pixel accuracy)
                u_c = sum(cols .* intensities) / total;
                v_c = sum(rows .* intensities) / total;

                stars(k, :) = [u_c, v_c, total];
                valid(k)    = true;
            end

            stars = stars(valid, :);
        end

        function [labels, n] = label_connected_components(obj, mask)
            % labeling the near-8 components
            %
            % algorithm: iterativa flood-fill w/ stack

            [H, W] = size(mask);
            labels = zeros(H, W);
            n = 0;

            for r = 1:H                 % defines the domain of the rows (between 1 and height of camera sensor)
                for c = 1:W             % defines the domain of the columns (between 1 and width of camera sensor)
                    if mask(r, c) && labels(r, c) == 0    % if mask(r,c) is true and the labels matrix are both empty (signed by zeros)
                        n = n + 1;                        % numeration of blob increases each time
                        
                        % Flood fill da (r, c)
                        stack = [r, c];
                        while ~isempty(stack)              % flood-fill process: 
                            pr = stack(end, 1);            % 1. reads the last row in pr,pc
                            pc = stack(end, 2);
                            stack(end, :) = [];            % 2. removes that last row (POP)

                            if pr < 1 || pr > H || pc < 1 || pc > W, continue; end   % pixels rows and columns must (obv) fit in the frame dimensions
                            if ~mask(pr, pc) || labels(pr, pc) ~= 0, continue; end   % if mask is true and labels is zero (if the pixel is on and has not been read)

                            labels(pr, pc) = n;            % this pixel now is part of the n-blob

                            % adds the near-8 pixels
                            stack(end+1, :) = [pr-1, pc-1]; % here we PUSH the 8-near pixels at the end of the stack matrix, and in the for and whil cycle it renews the stack vector and iterates the flood-fill all over the fram until everything is labeled
                            stack(end+1, :) = [pr-1, pc  ];
                            stack(end+1, :) = [pr-1, pc+1];
                            stack(end+1, :) = [pr  , pc-1];
                            stack(end+1, :) = [pr  , pc+1];
                            stack(end+1, :) = [pr+1, pc-1];
                            stack(end+1, :) = [pr+1, pc  ];
                            stack(end+1, :) = [pr+1, pc+1];
                        end
                    end
                end
            end
        end

        function visualize(obj, img, stars)
            % draw the image with the overlapped blobs (for debugging).
            figure;
            imshow(img, []); hold on;
            plot(stars(:,1), stars(:,2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
            title(sprintf('Detected %d stars', size(stars,1)));
        end
    end
end