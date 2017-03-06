function checkDir(path2outputs)

if ~exist(path2outputs, 'dir')
    mkdir(path2outputs);
end