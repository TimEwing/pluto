unset PLUTO_SPECTRUM
unset CHARON_SPECTRUM

python manager.py pluto NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_pluto_435.json --o output/pluto_435.json
python manager.py charon NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_charon_435.json --o output/charon_435.json
python manager.py pluto NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_pluto_555.json --o output/pluto_555.json
python manager.py charon NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_charon_555.json --o output/charon_555.json

export PLUTO_SPECTRUM='data/spectra/stis_solar.dat'
export CHARON_SPECTRUM='data/spectra/stis_solar.dat'

python manager.py pluto NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_pluto_435.json --o output/pluto_435_solar.json
python manager.py charon NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_charon_435.json --o output/charon_435_solar.json
python manager.py pluto NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_pluto_555.json --o output/pluto_555_solar.json
python manager.py charon NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_charon_555.json --o output/charon_555_solar.json

export PLUTO_SPECTRUM='data/spectra/pholus_spectrum.dat'
export CHARON_SPECTRUM='data/spectra/pholus_spectrum.dat'

python manager.py pluto NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_pluto_435.json --o output/pluto_435_pholus.json
python manager.py charon NH_RED NH_BLUE HST_F435W --fourier data/buie_fourier_charon_435.json --o output/charon_435_pholus.json
python manager.py pluto NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_pluto_555.json --o output/pluto_555_pholus.json
python manager.py charon NH_RED NH_BLUE HST_F555W --fourier data/buie_fourier_charon_555.json --o output/charon_555_pholus.json

unset PLUTO_SPECTRUM
unset CHARON_SPECTRUM
