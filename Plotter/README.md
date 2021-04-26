# Plotter

A tool which generates plots for your own analysis. By default, the xml file "plots.xml" will be loaded. All the configuration between the tags &lt;drawing&gt;(...)&lt;/drawing&gt; will automatically be loaded and will generate pdf files. Check the file "plots.xml" as an example.

## Parameters

   * &lt;enable&gt;: 1 (plot), 0 (do not plot)
   * &lt;signal_file&gt;: root signal file.
   * &lt;bkg_file&gt;: root background file. if not specified, it needs to be "none"). Signal file will be used also as the background.
   * &lt;variable_name&gt;: name of the branch to be plotted.
   * &lt;title_and_axis_name&gt;: as root parameter.
   * &lt;bin_size&gt;: as root parameter.
   * &lt;first_bin&gt;: as root parameter.
   * &lt;last_bin&gt;: as root parameter.
   * &lt;signal_cuts&gt;: cuts for the signal (combination of branches). Please, use [ASCII for html](http://www.asciitable.com/) for your symbols (>, <, &).
   * &lt;bkg_cuts&gt;: cuts for the background (combination of branches). Please, use [ASCII for html](http://www.asciitable.com/) for your symbols (>, <, &).
   * &lt;tagname&gt;: name of the pdf file (it does not to specify .pdf). 
 
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
