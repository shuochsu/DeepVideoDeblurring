# params
ALIGNMENT="_OF"
MODEL="model2_symskip_nngraph2_deeper"
METHOD=""
REALSET="data/testing_real_all_nostab"$ALIGNMENT
FN="1024_"$MODEL$ALIGNMENT$METHOD
FNREAL=$FN"_real"
GPUID=0
EPOCHTEST=400

export CUDA_VISIBLE_DEVICES=$GPUID

#prediction on qualitative (real) set
for subset in 'alley' 'anita' 'bicycle' 'boat' 'books' 'bridge' 'camp' 'car' 'cutting' 'kid' 'metro' 'office' 'opera' 'piano' 'playground' 'quinault1' 'starbucks' 'starbucks2' 'street' 'street2' 'sup' 'walk'
do
	th inference.lua -g $GPUID --model $MODEL --data_root $REALSET/$subset --model_param logs/$FN/param_epoch_$EPOCHTEST.t7 --bn_meanstd logs/$FN/bn_meanvar_epoch_$EPOCHTEST.t7 --saveDir outImg/$FNREAL/$subset --start_id 1 --n 100 --patchbypatch 1 --patch_size 256
done
