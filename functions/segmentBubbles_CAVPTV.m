function bubbleMask = segmentBubbles_CAVPTV(frameGray,parameters)
    
luminescence = parameters.xrayill;

if luminescence ==1
    % Bubble segmentation steps
    A = frameGray;
    A1 = medfilt2(A, [2 2]);
    edges = edge(A1, 'canny', [0.1 0.4]);
    contrastEnhanced = imadjust(uint8(edges) * 255);
    thresholded = contrastEnhanced > 0;
    se = strel('disk', 3);
    closedImage = imclose(thresholded, se);
    filledImage = imfill(closedImage, 'holes');
    areaThreshold = 50;
    filledImage = bwareaopen(filledImage, areaThreshold);
    se2 = strel('disk', 2);
    dilatedImage = imdilate(filledImage, se2);
    dilatedImage = imdilate(dilatedImage, se2);
    closedAgain = imclose(dilatedImage, se);
    filledAgain = imfill(closedAgain, 'holes');
    se3 = strel('disk', 2);
    erodedImage = imerode(filledAgain, se3);
    erodedImage = imerode(erodedImage, se3);
    finalMask = medfilt2(erodedImage, [10 10]);
    areaThreshold2 = 60;
    finalmask1 = bwareaopen(finalMask,areaThreshold2);
    bubbleMask = finalmask1;


elseif luminescence == 2
    % Bubble segmentation steps
    A = frameGray;
    A1 = medfilt2(A, [2 2]);
    edges = edge(A1, 'canny', [0.1 0.4]);
    contrastEnhanced = imadjust(uint8(edges) * 255);
    thresholded = contrastEnhanced > 0;
    se = strel('disk', 3);
    closedImage = imclose(thresholded, se);
    filledImage = imfill(closedImage, 'holes');
    areaThreshold = 50;
    filledImage = bwareaopen(filledImage, areaThreshold);
    se2 = strel('disk', 2);
    dilatedImage = imdilate(filledImage, se2);
    dilatedImage = imdilate(dilatedImage, se2);
    closedAgain = imclose(dilatedImage, se);
    filledAgain = imfill(closedAgain, 'holes');
    se3 = strel('disk', 2);
    erodedImage = imerode(filledAgain, se3);
    erodedImage = imerode(erodedImage, se3);
    finalMask = medfilt2(erodedImage, [10 10]);
    areaThreshold2 = 60;
    finalmask1 = bwareaopen(finalMask,areaThreshold2);
    bubbleMask = finalmask1;


elseif luminescence == 3 %when no particles and just inception (just for bubbles)
    A = frameGray;
    A1 = medfilt2(A, [2 2]);
    edges = edge(A1, 'canny', [0.1 0.4]);
    contrastEnhanced = imadjust(uint8(edges) * 255);
    thresholded = contrastEnhanced > 0;
    se = strel('disk', 3);
    closedImage = imclose(thresholded, se);
    filledImage = imfill(closedImage, 'holes');
    areaThreshold = 50;
    filledImage = bwareaopen(filledImage, areaThreshold);
    se3 = strel('disk', 1);
    erodedImage = imerode(filledImage, se3);
    bubbleMask = erodedImage;
    
elseif luminescence == 4 %when no particles and just inception [select all size of bubbles]
    A = frameGray;
    % A1 = medfilt2(A, [2 2]);
    % figure
    % imshow(A1); title('med filt');

    % Lower the thresholds for edge detection to be more sensitive
     edges = edge(A, 'canny', [0.05 0.2]);

    % Contrast enhancement (you might also consider working directly with 'edges' if desired)
    contrastEnhanced = imadjust(uint8(edges) * 255);
    thresholded = contrastEnhanced > 0;

    % Use a smaller structuring element for closing to better preserve small features
    se = strel('disk', 3);
    closedImage = imclose(thresholded, se);

    % Fill any holes inside the bubbles
    filledImage = imfill(closedImage, 'holes');

    % Lower the area threshold so that small bubbles are not removed
    areaThreshold = 2;  % adjust as needed based on your bubble sizes
    filledImage = bwareaopen(filledImage, areaThreshold);

    % Optionally, reduce or skip erosion to avoid shrinking small bubbles
    % Here, either skip erosion or use a very small structuring element.
    % For example, if you decide erosion is still needed:
    se3 = strel('disk', 1);  % effectively no erosion; you can try disk(1) if needed
    erodedImage = imerode(filledImage, se3);

    bubbleMask = erodedImage;

    % Fill any holes inside the bubbles
    bubbleMask = imfill(bubbleMask, 'holes');

    %open 
    se4 = strel('disk',1);
    bubbleMask = imopen(bubbleMask,se4);


elseif  luminescence == 5 % Cloud cavitation case
    A= frameGray;
    bubbleMask = adaptthresh(A);


else
    error('Unsupported luminescence value: %d, please modify the parameters.xrayill value', luminescence);
end


end