python main.py pluto NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_pluto_435.json --o output/pluto_435.json
python main.py pluto NH_RED NH_BLUE HST_F435W --errorbar upper --fourier data/buie_fourier_pluto_435.json --o output/pluto_435_upper.json
python main.py pluto NH_RED NH_BLUE HST_F435W --errorbar lower --fourier data/buie_fourier_pluto_435.json --o output/pluto_435_lower.json

python main.py charon NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_charon_435.json --o output/charon_435.json
python main.py charon NH_RED NH_BLUE HST_F435W --errorbar upper --fourier data/buie_fourier_charon_435.json --o output/charon_435_upper.json
python main.py charon NH_RED NH_BLUE HST_F435W --errorbar lower --fourier data/buie_fourier_charon_435.json --o output/charon_435_lower.json


python main.py pluto NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_pluto_555.json --o output/pluto_555.json
python main.py pluto NH_RED NH_BLUE HST_F555W --errorbar upper --fourier data/buie_fourier_pluto_555.json --o output/pluto_555_upper.json
python main.py pluto NH_RED NH_BLUE HST_F555W --errorbar lower --fourier data/buie_fourier_pluto_555.json --o output/pluto_555_lower.json

python main.py charon NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_charon_555.json --o output/charon_555.json
python main.py charon NH_RED NH_BLUE HST_F555W --errorbar upper --fourier data/buie_fourier_charon_555.json --o output/charon_555_upper.json
python main.py charon NH_RED NH_BLUE HST_F555W --errorbar lower --fourier data/buie_fourier_charon_555.json --o output/charon_555_lower.json
