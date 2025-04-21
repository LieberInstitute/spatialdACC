 fname = '/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/if-images/V12N28-333_D1.mat';
 img = load(fname,'DAPI');
 imwrite(mat2gray(img.DAPI(5000:8000, 4000:7000)),'/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/tile1_D1.png');
 imwrite(mat2gray(img.DAPI),'/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/D1.png');
 
 fname = '/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/if-images/V12N28-333_A1.mat';
 img = load(fname,'DAPI');
 imwrite(mat2gray(img.DAPI(8000:11000, 4000:7000)),'/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/tile1_A1.png');
 imwrite(mat2gray(img.DAPI),'/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/A1.png');
 