# ATLAS MC -> HyPER training inputs

This repository holds the code to convert ATLAS datasets ( produced with TCT ) into .h5 files. The code uses a FastFrames custom class that supports various channels, at the moment the supported channels are:
- TtbarAllHadronic
- TtbarSingleLepton

After FastFrames produces the desired ROOT files one needs to run the corresponding ROOT -> H5 converter under `H5Converters`.