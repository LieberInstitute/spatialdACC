%copyfile /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/if-images/*1.tif /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/samui_input/

% attach segmented image to raw image
myfiles = dir('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/cellpose_masks/*1_cp_masks.png');
dr = '/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/samui_input/';
for i = 1:numel(myfiles)
fname=fullfile(myfiles(i).folder, myfiles(i).name);
BW = imread(fname);
BW = label2rgb(BW,'hsv','w','shuffle');
imwrite(BW,fullfile(dr, [myfiles(i).name(1-19),'.tif']),"WriteMode","append")
disp(i)
end

