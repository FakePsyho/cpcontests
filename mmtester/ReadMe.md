Simple local tester for Topcoder's Marathon contests. 

# Setup
- Install tabulate: `pip install tabulate`
- Copy `tester.py` and `tester.cfg` to the directory with the problem
- Edit `tester.cfg` in order to set the number of threads, tester jar file, executable name & scoring function

# Main Features
- Runs batches of tests and stores results in a JSON format
- Grouping/Filtering results by custom metadata
- Searching through tests according to some metrics (useful for finding tests close to timing out or some more extreme test cases)

# How it Works
- Run `tester.py [results_file]` in order to run a batch of tests and create a `results_file.res` with scores. Output of each test is redirected to a text file and then later the file is scanned in order to extract the score (and potentially other variables). Results file will contain scores for each test case in JSON format. By default, output files will be stored at `./tests/results_file`.
- Run `tester.py -s` to see the scoreboard for all of the results' files in the current directory. 
- By default, results files contains only Score and the test ID (seed). You can add additional metadata by printing `[DATA] X = Y` (spaces are mandatory) in your solution. For example, let's say the the problem is on the square grid of size N. If your solution prints `[DATA] N = ...`, you can run `tester.py -s --data [results_file_with_data] --groups N` to see the summary grouped by different sizes of the grid. You can also manually specify multiple groups with different criteria like `--groups N=2 N=3-5 N=6-10 D=2 D=4`. Similarly, you can use `--filters N=5-10` to filter the testscases to include only those fullfilling specific criteria.
- For more options/features, please look at the config comments and/or built-in help. You can override most of the options by specifying adequate command line option.

# Notes
- There's a frequent problem with marathon testers where stderr (debug stuff) and stdout (official output from the tester) is mangled together. If it happens then `tester.py` will be unable to read the score and unfortunately the whole run will fail. There's currently no nice workaround for this and the best way to avoid this problem is to edit the source code of the java tester and either add locking mechanism on output or add some sleep just before printing the score.
- Feel free to report any errors and/or feature requests (DMing me on twitter might be the easiest way of getting my attention). Can't promise anything, but I might consider fixing/adding those.
- It should be easy to extend to other contests with similar structure (i.e. single-player problems). You will have to edit the `run_test` function.
- Tested with python 3.7.6 on Windows & Linux.
