% For more information, please refer to the following papers:
% Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% This code is written by Mingxing Zhang, University of Virginia.
function write3Dtiff_V2(array, filename)

outtiff=Tiff(filename, 'w');

dims = size(array);
if length(dims)==2
    dims(3)=1;
end

tagstruct.ImageLength = dims(1);
tagstruct.ImageWidth = dims(2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if isa(array, 'single')
    tagstruct.BitsPerSample = 32;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
elseif isa(array, 'uint16')
    tagstruct.BitsPerSample = 16;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
else
    tagstruct.BitsPerSample = 8;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
end
tagstruct.SamplesPerPixel = 1;
%tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

%outtiff.setTag(tagstruct)

for z=1:dims(3)
    %outtiff.currentDirectory
    outtiff.setTag(tagstruct)
%     outtiff.setTag('ImageWidth', dims(1))
%     outtiff.setTag('ImageLength', dims(2))
%     outtiff.setTag('Photometric',rawtiff.getTag('Photometric'))
%     outtiff.setTag('PlanarConfiguration',rawtiff.getTag('PlanarConfiguration'))
%     outtiff.setTag('BitsPerSample',32)
%     outtiff.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
    outtiff.write(array(:,:,z))
    outtiff.writeDirectory()
end
outtiff.close()
