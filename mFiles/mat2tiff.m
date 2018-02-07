function mat2tiff(vol,saveslice,numclass, filename)
% mat2tiff(vol,saveslice,numclass,filename)
% mat2tiff converts a .mat file to a .tiff (uint16).  This works for 2D (image) data
%   and 3D (image stack/volume) data.
%
%
%INPUTS
%------
% vol: volume you want to save
% saveslice: in what slice direction do you want to save the data
%               (options 'xy','xz',& 'yz')
% numclass: numeric class of output data (options: 'uint8', 'uint16',
%                                                 'double', & 'same')
% filename: filename for the volume you want to save (default will open
% prompt).
%
% Author: Eyal Bar-Kochba, Brown University (2011) <><

%
if exist('filename','var') == 0
    prompt = {'Enter the filename for the tiff stack',};
    filename = inputdlg(prompt,'',1,{'name'});
    filename = cell2mat(filename);
end

filename = [filename,'.tif'];

switch lower(saveslice)
    case 'xz', vol = permute(vol, [2 3 1]);
    case 'yz', vol = permute(vol, [1 3 2]);
    case 'xy', vol = permute(vol, [1 2 3]); 
end

switch lower(numclass)
    case 'uint16', vol = uint16(vol);
    case 'uint8',  vol = im2uint8(vol);
    case 'double', vol = double(vol);
end


if exist(filename,'file')
    choice = questdlg(['Would you like to overwrite', num2str(filename),'?'], ...
        'File exists in current directory', ...
        'Yes','No','Cancel','Cancel');
    switch choice
        case 'Yes',
            delete(filename); 
            save3D(vol,filename,saveslice)
        case 'No',
            prompt = {'Enter the filename for the tiff stack',};
            filename = inputdlg(prompt,'',1,{'name'});
            filename = cell2mat(filename);
            save3D(vol,filename,saveslice)
            
        otherwise,
    end
else
    
    saveName = fullfile(pwd,filename);
    save3D(vol,saveName,saveslice)
end

function save3D(vol,saveName,saveslice)
wb = waitbar(0,'');
for ind = 1:size(vol,3)
    waitbar(ind/size(vol,3),wb,['Saving ',saveslice,': ',num2str(ind),'/',num2str(size(vol,3))])
    imwrite(vol(:,:,ind),saveName,'WriteMode','append','Compression','none');
end
close(wb);