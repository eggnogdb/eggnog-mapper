/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --master --cport 51400 --wport 51401 --hmmdb hmmdb/euk_500/euk_500.hmm --daemon --pid ./euk.m.pid
/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --master --cport 51500 --wport 51501 --hmmdb hmmdb/bact_50/bact_50.hmm --daemon --pid ./bact.m.pid
/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --master --cport 51600 --wport 51601 --hmmdb hmmdb/arch_1/arch_1.hmm   --daemon --pid ./arch.m.pid

/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --worker --wport 51401 --hmmdb hmmdb/euk_500/euk_500.hmm --daemon --pid ./euk.w.pid  --cpu 20
/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --worker --wport 51501 --hmmdb hmmdb/bact_50/bact_50.hmm --daemon --pid ./bact.w.pid --cpu 20
/kappa/data/eggnog_deps/bin/hmmer-3.1b2/src/hmmpgmd --worker --wport 51601 --hmmdb hmmdb/arch_1/arch_1.hmm   --daemon --pid ./arch.w.pid --cpu 20
