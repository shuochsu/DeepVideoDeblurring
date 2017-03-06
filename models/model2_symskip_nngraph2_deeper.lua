require 'nn'
-- require 'cunn'
require 'nngraph'

local Convolution = nn.SpatialConvolution
local UpConvolution = nn.SpatialFullConvolution
local SBatchNorm = nn.SpatialBatchNormalization
local Max = nn.SpatialMaxPooling
local Avg = nn.SpatialAveragePooling
local ReLU = nn.ReLU

function create_model()
	input = - nn.Identity()
	ref = input 
		 - nn.Narrow(2, 7, 3)
	F0 =  input 
		 - Convolution(15,64,5,5,1,1,2,2)
	     - SBatchNorm(64,1e-3)
	     - ReLU(true)

    D1 =  F0 
		 - Convolution(64,64,3,3,2,2,1,1)
	     - SBatchNorm(64,1e-3)
	     - ReLU(true)
	F1 =  D1 
		 - Convolution(64,128,3,3,1,1,1,1)
	     - SBatchNorm(128,1e-3)
	     - ReLU(true)
	F2 =  F1 
		 - Convolution(128,128,3,3,1,1,1,1)
	     - SBatchNorm(128,1e-3)
	     - ReLU(true)

	D2 =  F2 
		 - Convolution(128,256,3,3,2,2,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)
	F3 =  D2
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)
	F4 =  F3
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)
	F5 =  F4
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)

	D3 =  F5 
		 - Convolution(256,512,3,3,2,2,1,1)
	     - SBatchNorm(512,1e-3)
	     - ReLU(true)
	F6 =  D3
		 - Convolution(512,512,3,3,1,1,1,1)
	     - SBatchNorm(512,1e-3)
	     - ReLU(true)
	F7 =  F6
		 - Convolution(512,512,3,3,1,1,1,1)
	     - SBatchNorm(512,1e-3)
	     - ReLU(true)
	F8 =  F7
		 - Convolution(512,512,3,3,1,1,1,1)
	     - SBatchNorm(512,1e-3)
	     - ReLU(true)

	U1 =  F8 
		 - UpConvolution(512,256,4,4,2,2,1,1)
	     - SBatchNorm(256,1e-3)
    S1 = {F5,U1}
		 - nn.CAddTable()
	 	 - ReLU(true)
	F9 =  S1
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)
	F10 =  F9
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)
	F11 =  F10
		 - Convolution(256,256,3,3,1,1,1,1)
	     - SBatchNorm(256,1e-3)
	     - ReLU(true)

	U2 =  F11
		 - UpConvolution(256,128,4,4,2,2,1,1)
	     - SBatchNorm(128,1e-3)
	S2 = {F2,U2}
		 - nn.CAddTable()
	 	 - ReLU(true)
	F12 =  S2 
		 - Convolution(128,128,3,3,1,1,1,1)
	     - SBatchNorm(128,1e-3)
	     - ReLU(true)
	F13 =  F12 
		 - Convolution(128,64,3,3,1,1,1,1)
	     - SBatchNorm(64,1e-3)
	     - ReLU(true)

	U3 =  F13 
		 - UpConvolution(64,64,4,4,2,2,1,1)
	     - SBatchNorm(64,1e-3)
    S3 = {F0,U3}
		 - nn.CAddTable()
	 	 - ReLU(true)
	F14 =  S3 
		 - Convolution(64,15,3,3,1,1,1,1)
	     - SBatchNorm(15,1e-3)
	     - ReLU(true)
	F15 =  F14 
		 - Convolution(15,3,3,3,1,1,1,1)
	     - SBatchNorm(3,1e-3)

	S4 = {ref,F15}
		 - nn.CAddTable()
	 	 - nn.Sigmoid()

	g = nn.gModule({input},{S4})

-- 	indata = torch.rand(4,15,64,64)
-- 	g:forward(indata)

-- 	graph.dot(g.fg, 'model2_symskip_deeper', 'model2_symskip_deeper')

	return g
end