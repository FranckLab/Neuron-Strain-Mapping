# Neuron Strain Mapping
The included code is a companion to the manuscript (In submission) and provides a method for transformation and resolution of both, locally discrete, and mean cellular deformation strains across single cells or networks given a known far-field strain. The code also includes an image-processing workflow, part of which takes advantage of the freely available NeuronStudio software to make neuron tree structures. 

#### System Recommendations:
* Intel Core i7-4790K CPU, 
* 16 GB memory, 
* Matlab version R2015b or newer,
* Windows 10

#### Required non-standard software:
* Matlab
* NeuronStudio (available at http://research.mssm.edu/cnic/tools-ns.html)

#### Optional non-standard software:
* TecPlot
* (other 3D viewers or Matlab's native viewer may be used.)

#### Installation guide:
Please refer to NeuronStudio documentation for related installation instructions.
The rest of the code does not need installing. Place it in your working directory.

#### Demo:
[Example Data](https://drive.google.com/drive/folders/1Em01Jwa_PbfkEjtRQMGAua3c6PMrK5dy?usp=sharing)
The data is already in .mat format. Since there is no needed bioformat to mat conversion, 
block 2 of ProcessRunFile can be skipped.
We recommend keeping the data in the folder matching in name because the code expects that architecture.
 
The expected output should be a processed .tiff file, an .swc file output from NeuronStudio, a segmented.tiff file,
and a .plt file for optional import into TecPlot

#### Expected Run Time:
Expected run time varies by dataset. Processing of included demo data set on 
a computer with recommended system requirements resulted in the following run times:

* ProcessRunFile: 0.7s (without FIDVC drift correction)
* MaskRunFile: 4.1s 
* StrainRunFile: 0.9s (without plotting)
* BlebStrainRunFile: 1.6s

The FIDVC (optional block, code not included. Download FIDVC code from our [Github page](https://github.com/FranckLab/FIDVC)) 
and bfmatlab folders need to be in the search path in order to run all blocks of ProcessFunFile. 

#### Contact and support
Please refer to the manuscript text or comments within the code for detailed instructions on how to run it. Also refer to [Questions/Issues](https://github.com/FranckLab/Neuron-Strain-Mapping/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).




