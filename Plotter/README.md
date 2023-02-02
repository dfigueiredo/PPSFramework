# Plotter

A tool which generates plots for your own analysis. By default, the xml file "plots.xml" will be loaded. All the configuration between the tags &lt;drawing&gt;(...)&lt;/drawing&gt; will automatically be loaded and will generate pdf files. Check the file "plots.xml" as an example.

## XML Parameters

| Options       | Comments |
| ------------- | -------------:|
| enable | 1 (plot), 0 (do not plot) |
| signal_file | root signal file |
| bkg_file | root background file. If not needed, it must be set  "none" |
| variable_name | name of the branch to be plotted |
| title_and_axis_name | as root parameter |
| bin_size | as root histogram plotter parameter |
| first_bin | as root histogram plotter parameter |
| last_bin | as root histogram plotter parameter |
| signal_cuts | cuts for the signal (combination of branches). Please, use [ASCII for html](http://www.asciitable.com/) for your symbols (>, <, &) |
| bkg_cuts | cuts for the background (combination of branches). Please, use [ASCII for html](http://www.asciitable.com/) for your symbols (>, <, &) |
| tagname | name of the output pdf file (no need to specify .pdf) |

As an example, check the file https://github.com/dfigueiredo/PPSFramework/blob/main/Plotter/plots.xml

```sh
python PlotterTool.py
draw --file your_file.xml [press enter]
```

For more options:

```sh
python PlotterTool.py --h
Usage: PlotterTool.py [options]

Options:
  -h, --help            show this help message and exit
  -f FILE, --filename=FILE
                        XML mapping file
  -p PARSING, --parsing=PARSING
                        parsing: commands which can be passed from SHELL
                        directly. [parsing: --p "draw --file filename.xml"]
  -v, --verbose         make lots of noise [default]
```

**Important**: all the histograms cosmetics should be changed in the Plotter.cc source code. Remember to compile it, otherwise the python frontend interface will not work.

```sh
make clean
make
```

# Event List

A tool to generate a file with the selected events after a pre-selected cut.

```sh
make clean
make
./EventList --f tree_pps.root
```
