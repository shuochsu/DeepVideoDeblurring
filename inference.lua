require 'xlua'
require 'optim'
require 'cunn'
require 'image'
require 'gnuplot'

local c = require 'trepl.colorize'

opt = lapp[[
    -g, --gpuid             (default 0)                  gpu id
    --seed                  (default 123)                rand seed

    --model                 (default 'model')            model name
    --model_param           (default 'param')            weight file
    --bn_meanstd            (default 'bn_meanvar')       batch normalization file

    --num_frames            (default 5)                  number of frames in input stack
    --num_channels          (default 3)                  rgb input
    --max_intensity         (default 255)                maximum intensity of input
    --patchbypatch          (default 0)                  tile-based predication
    --patch_size            (default 128)                should be equal to that during training

    --data_root             (default 'data/testing')     folder for testing data
    --saveDir               (default 'debug')            folder for intermediate prediction result
    
    --start_id              (default 1)                  start from which file id
    --n                     (default 1)                  predict on how many frames

]]

local function load_model()
    require('models/' .. opt.model .. '.lua')

    pred = create_model(opt.num_frames*opt.num_channels):cuda()

    local param, grad_param = pred:getParameters()
    print(string.format('number of parameters: %d', param:nElement()))
    
    print(c.blue '==>' ..' loading parameters')
    -- load parameters
    local params = torch.load(opt.model_param)
    assert(params:nElement() == param:nElement(), string.format('%s: %d vs %d', 'loading parameters: dimension mismatch.', params:nElement(), param:nElement()))
    param:copy(params)

    if(string.len(opt.bn_meanstd) > 0) then 
        local bn_mean, bn_std = table.unpack(torch.load(opt.bn_meanstd))

        for k,v in pairs(pred:findModules('nn.SpatialBatchNormalization')) do
            v.running_mean:copy(bn_mean[k])
            v.running_var:copy(bn_std[k])
        end
        pred:evaluate() -- bn statistics required
    end
end

function my_forward( m, img )    
    local data = img:view(1, img:size(1), img:size(2), img:size(3))
    out = m:forward(data)
    return out
    -- for i = 1, #m do
    --     data = m.modules[i]:updateOutput(data)
    --     -- print('..' .. i)
    --     if m.modules[i].finput then
    --         m.modules[i].finput:set()
    --     end
    -- end
    -- return data:clone()
end

print(opt)

torch.manualSeed(opt.seed)

lleft_folder = 'image_-2'
left_folder = 'image_-1'
mid_folder = 'image_0'
right_folder = 'image_1'
rright_folder = 'image_2'

paths.mkdir(opt.saveDir)
filename = paths.concat(opt.saveDir, string.format('opt.t7'))
torch.save(filename, opt)

-- load model
print(c.blue '==>' ..' configuring model')
load_model()

for i = opt.start_id, opt.start_id+opt.n-1 do
    local file_id = i
    print('--- ' ..  i .. ' --- ' .. ' fn: ' .. file_id)

    local global_mean = 0
    local global_std = 1

    local ll_fn = string.format('%s/%s/%05d.jpg', opt.data_root, lleft_folder, file_id)
    local l_fn = string.format('%s/%s/%05d.jpg', opt.data_root, left_folder, file_id)
    local m_fn = string.format('%s/%s/%05d.jpg', opt.data_root, mid_folder, file_id)
    local r_fn = string.format('%s/%s/%05d.jpg', opt.data_root, right_folder, file_id)
    local rr_fn = string.format('%s/%s/%05d.jpg', opt.data_root, rright_folder, file_id)
    local ll_img = image.load(ll_fn, opt.num_channels, 'byte'):float()
    local l_img = image.load(l_fn, opt.num_channels, 'byte'):float()
    local m_img = image.load(m_fn, opt.num_channels, 'byte'):float()
    local r_img = image.load(r_fn, opt.num_channels, 'byte'):float()
    local rr_img = image.load(rr_fn, opt.num_channels, 'byte'):float()

    ll_img:div(opt.max_intensity):add(-global_mean):div(global_std)
    l_img:div(opt.max_intensity):add(-global_mean):div(global_std)
    m_img:div(opt.max_intensity):add(-global_mean):div(global_std)
    r_img:div(opt.max_intensity):add(-global_mean):div(global_std)
    rr_img:div(opt.max_intensity):add(-global_mean):div(global_std)

    local img_c = l_img:size(1)
    local img_h = l_img:size(2)
    local img_w = l_img:size(3)

    local stack = torch.cat({ll_img, l_img, m_img, r_img, rr_img}, 1):cuda()

    if opt.patchbypatch == 0 then
        -- pad if not divideble by 8
        local hP = 8*torch.ceil(img_h/8)
        local wP = 8*torch.ceil(img_w/8)
        local stackP = torch.Tensor(stack:size(1),hP,wP):zero()
        stackP[{{},{1,img_h},{1,img_w}}] = stack:double()
        stackP = stackP:cuda()
        out = my_forward(pred, stackP)
        out = out:view(opt.num_channels, hP, wP)
        out = out[{{},{1,img_h},{1,img_w}}]
    else
        signal = require 'signal'
        local patch_shift = torch.floor(opt.patch_size/2)
        local hP = patch_shift*torch.ceil((img_h-opt.patch_size)/patch_shift) + opt.patch_size
        local wP = patch_shift*torch.ceil((img_w-opt.patch_size)/patch_shift) + opt.patch_size
        local stackP = torch.Tensor(stack:size(1),hP,wP):zero()
        stackP[{{},{1,img_h},{1,img_w}}] = stack:double()
        stackP = stackP
        local image_deblurred = torch.Tensor(img_c,hP,wP):zero()
        local image_norm = torch.Tensor(img_c,hP,wP):zero()
        local window = signal.hann(opt.patch_size):resize(1,opt.patch_size)
        local window_t = window:t()
        window = torch.mm(window_t,window)
        window = torch.add(torch.repeatTensor(window,img_c,1,1), torch.Tensor(img_c,opt.patch_size,opt.patch_size):fill(1e-6))

        for j = 1, hP-opt.patch_size+1, patch_shift do
            for k = 1, wP-opt.patch_size+1, patch_shift do
                
                local jmax =  torch.min(torch.Tensor({j+opt.patch_size-1,hP}));
                local kmax =  torch.min(torch.Tensor({k+opt.patch_size-1,wP}));
                
                local stack_crop = stackP[{{}, {j,jmax}, {k,kmax}}]:contiguous()

                local curr_patch_deblurred = my_forward(pred, stack_crop:cuda())
                curr_patch_deblurred = curr_patch_deblurred:view(opt.num_channels, opt.patch_size, opt.patch_size)

                image_deblurred[{{}, {j,jmax}, {k,kmax}}] = image_deblurred[{{}, {j,jmax}, {k,kmax}}] + torch.cmul(curr_patch_deblurred:double(), window);
                image_norm[{{}, {j,jmax}, {k,kmax}}] = image_norm[{{}, {j,jmax}, {k,kmax}}] + window;
            end
        end

        image_deblurred = torch.cdiv(image_deblurred, image_norm);
        out = image_deblurred[{{}, {1,img_h}, {1,img_w}}]
    end

    out:mul(global_std):add(global_mean):mul(opt.max_intensity)
    image.save(string.format('%s/%05d.jpg', opt.saveDir, file_id), out:byte())

    print('writing deblurred image done..')

end





