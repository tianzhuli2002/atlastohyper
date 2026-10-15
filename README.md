# ATLAS MC -> HyPER training inputs

This repository holds the code to convert ATLAS datasets ( produced with TCT ) into .h5 files. The code uses a FastFrames custom class that supports various channels, at the moment the supported channels are:
- TtbarAllHadronic
- TtbarSingleLepton

After FastFrames produces the desired ROOT files one needs to run the corresponding ROOT -> H5 converter under `H5Converters`.

## Once you have H5 files

Use the [HyPER repository](https://github.com/tzuhanchang/HyPER) to train your models.

## Once your models are trained

Use the following repository to [export your pytorch models to ONNX](https://github.com/diegobaronm/HyPER_ONNX).

Congratulations! You can use your `.onnx` files to perform inference in [TopCPToolkit](https://topcptoolkit.docs.cern.ch) or [FastFrames](https://atlas-project-topreconstruction.web.cern.ch/fastframesdocumentation/latest/).# atlastohyper
