# experiment name
exp=disc_cov_dep

# sample size and signal strength
n1=50
n2=50
lambda=15

# number of trials
ntrials=50

# p=51
p=51
now=$(date '+%Y-%m-%d_%H:%M')
out_name=$exp\_ntrials$ntrials\_p$p\_n1\_$n1\_n2\_$n2\_lambda$lambda\_$now.Rout
printf "\n\nRunning $out_name...\n"
SECONDS=0
R CMD BATCH --no-save --no-restore "--args \
experiment='$exp' \
p=$p \
n1=$n1 \
n2=$n2 \
lambda=$lambda \
n_trials=$ntrials" \
skips=c('mgm')"
test.R $out_name
duration=$SECONDS
echo "Done. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"

# p=32

# p=11
