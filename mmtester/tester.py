#!/usr/bin/env python3

# Author: Psyho
# Twitter: https://twitter.com/fakepsyho

#TODO:
# HIGH PRIORITY:
# -check for config in current dir otherwise use default & add tool to bin
# -add parameter to use tests/run_name/. instead of tests/.
# -fix double printing progress (HOW?)
# -fix problem when missing score (HOW?)
# -add auto buckets for groups like D(3) or D=2-10(4)
# LOW PRIORITY:
# -use --tests for --find?
# -add wrapper for fatal errors
# -show: print the # of tests for each group (line below the header?)
# -show: add transpose
# -more error checking / clearer error messages
# -add comments to code (explaining functions should be enough)
# -add more annotations to functions
# -add "--data LATEST" so that newest res file is always used for metadata
# -fix grouping/filtering if data file doesn't contain all test cases
# -add warnings if data is not present for all test cases?
# -add support for custom scoring (cfg would be python code?)
# -add RUNNER parameters (like for hash code) (Moved to RUNNER?)
# -add batching? (Moved to RUNNER?)
# -sync with RUNNER? (how?)
# -add cleanup on ctrl+c (what that would be?)
# -add to pypi
# -automatic way to install dependencies (pip install tabulate + ??)
# -change to subparsers (exec / show / find?)
# ???:
# -is cfg[...][...] really the best way to store config? converting types is really ugly :(
# -is it possible to shrink the column name with tabulate (cur min is header width + 2 spaces)


import tabulate
import re
import math
import sys
import os
import argparse
import subprocess
import glob
import json
import time
import configparser
from typing import List, Dict, Union
import queue
from threading import Thread

global args
args = None

global cfg
cfg = None

DEFAULT_CONFIG_PATH = 'tester.cfg'

tests_queue = queue.Queue()
results_queue = queue.Queue()


def try_str_to_numeric(x):
    if x is None:
        return None

    try: 
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


def run_test(seed) -> Dict:
    run_dir = args.name or '_default'
    output_path = f'{cfg["general"]["tests_dir"]}/{run_dir}/seed{seed}.out'
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    
    subprocess.run(f'{cfg["general"]["run_cmd"]} -exec "{args.exec}" -seed {seed} {args.tester_arguments} > {output_path}', shell=True)
    rv = {'id': seed}
    with open(output_path) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) == 3 and tokens[0] == 'Score' and tokens[1] == '=':
                rv['score'] = float(tokens[2])
            if len(tokens) == 4 and tokens[0] == '[DATA]' and tokens[2] == '=':
                rv[tokens[1]] = try_str_to_numeric(tokens[3])
                
    if 'score' not in rv:
        print(f'\r[Error] Seed: {seed} cointains no score')
        
    if cfg['general']['keep_output_files'].lower() != 'true':
        os.remove(output_path)
    return rv
    
    
def worker():
    while True:
        try:
            seed = tests_queue.get(False)
        except queue.Empty:
            return
        result = run_test(seed)
        results_queue.put(result)
    
    
def find_res_files(dir='.'):
    return glob.glob(f'{dir}/*{cfg["general"]["results_ext"]}')
    
    
def load_res_file(path) -> Dict[int, float]:
    if not os.path.exists(path):
        print(f'Fatal Error: Cannot locate {path} results file')
        sys.exit(1)

    with open(path) as f:
        lines = f.read().splitlines()
    results = [json.loads(line) for line in lines]
    return {result['id']: result for result in results} 


def process_raw_scores(scores: List[float], scoring: str): 
    if scoring=='raw':
        return scores
    if scoring=='min':
        best_score = min([math.inf] + [score for score in scores if score >= 0])
        return [best_score / score if score >= 0 else 0 for score in scores]
    if scoring=='max':
        best_score = max([0] + [score for score in scores if score >= 0])
        return [score / best_score if score > 0 else 0 for score in scores]
    
    print(f'Fatal Error: Unknown scoring function: {scoring}')
    sys.exit(1)
        
    
def apply_filter(tests, data, filter):
    var, range = filter.split('=')
    if '-' in range:
        lo, hi = range.split('-')
        lo = try_str_to_numeric(lo)
        hi = try_str_to_numeric(hi)
        return [test for test in tests if lo <= data[test][var] <= hi]
    else:
        value = try_str_to_numeric(range)
        return [test for test in tests if data[test][var] == value]

    
def show_summary(runs: Dict[str, Dict[int, float]], tests: Union[None, List[int]] = None, data=None, groups=None, filters=None):
    if not tests:
        tests_used = [set(run_results.keys()) for run_name, run_results in runs.items()]
        tests = tests_used[0].intersection(*tests_used[1:])
    else:
        # TODO: error check if tests are cointained in intersection of all results files?
        pass

    if not data and (filter or groups):
        print('[Error] Filter used but no data is provided')
        sys.exit(1)
        
    if filters:
        for filter in filters:
            tests = apply_filter(tests, filter)
            
    group_names = []
    group_tests = []
    
    group_names.append('Score')
    group_tests.append(tests)
    if groups:
        for group in groups:
            if '=' in group:
                group_names.append(group)
                group_tests.append(apply_filter(tests, data, group))
            else:
                var = group
                var_set = sorted(set([data[test][var] for test in tests]))
                for value in var_set:
                    group_names.append(f'{var}={value}')
                    group_tests.append(apply_filter(tests, data, f'{var}={value}'))
        
    headers = ['Run'] + group_names + ['Fails']
        
    table = [[run_name] for run_name in runs]
    
    total_fails = {run_name: 0 for run_name in runs}
    for group_test in group_tests:
        total_scores = {run_name: 0.0 for run_name in runs}
        for test in group_test:
            scores = process_raw_scores([run_results[test]['score'] for run_results in runs.values()], args.scoring)
            for run_name, score in zip(runs.keys(), scores):
                total_scores[run_name] += score
                total_fails[run_name] += 1 if score <= 0 else 0
        for i, run_name in enumerate(runs):
            table[i].append(total_scores[run_name] * (args.scale / len(group_test) if args.scale else 1.0))
    
    for i, run_name in enumerate(runs):
        table[i].append(total_fails[run_name])
    
    print('Total Tests:', len(tests))
    
    if args.scale:
        total_scores = {run_name: score * args.scale / len(tests) for run_name, score in total_scores.items()}
    longest_name = max([len(run_name) for run_name in runs])
        
    print(tabulate.tabulate(table, headers=headers))
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tester for Marathon Matches')
    parser.add_argument('name', type=str, nargs='?', default=None, help='name of the run') 
    parser.add_argument('-c', '--config', type=str, default=DEFAULT_CONFIG_PATH, help='path to cfg file')
    parser.add_argument('-t', '--tests', type=str, help='number of tests to run, range of seeds (e.g. A-B) or the name of the JSON/text file with the list of seeds')
    parser.add_argument('-m', '--threads_no', type=int, help='number of threads to use') 
    parser.add_argument('-e', '--exec', type=str, default=None, help='executable for the tester') 
    parser.add_argument('-p', '--progress', action='store_true', help='shows current progress when testing') 
    parser.add_argument('-a', '--tester_arguments', type=str, default='', help='additional arguments for the tester')
    parser.add_argument('-b', '--benchmark', type=str, default=None, help='benchmark res file to test against')
    parser.add_argument('-s', '--show', action='store_true', help='shows current results') 
    parser.add_argument('--data', type=str, default=None, help='file with metadata, used for grouping and filtering; in order to always use latest results file set it to LATEST') 
    parser.add_argument('--filters', type=str, default=None, nargs='+', help='filters results based on criteria') 
    parser.add_argument('--groups', type=str, default=None, nargs='+', help='groups results into different groups based on criteria') 
    parser.add_argument('--scale', type=float, help='sets scaling of results') 
    parser.add_argument('--scoring', type=str, default=None, help='sets the scoring function used for calculating ranking')
    parser.add_argument('--sorting', type=str, default=None, choices=['name', 'date'], help='sets how the show runs are sorted')
    parser.add_argument('--find', type=str, default=None, nargs='+', help='usage: --find res_file var[+/-] [limit]; sorts tests by var (asceding / descending) and prints seeds; can be combined with --filters')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.config):
        print(f'Fatal Error: Cannot locate config file at {args.config}. Check XXX for more information') #TODO: add a link to github page later / add option to automatically generate config file?
        sys.exit(1)
        
    cfg = configparser.ConfigParser()
    cfg.read(args.config)
    
    #XXX: probably there's a better to do this
    def convert(value, type=str):
        if value is None or value == '':
            return None
        
        if type == bool:
            return value.lower() in ['true', 'yes']
        return type(value)
        
    args.tests = try_str_to_numeric(args.tests or convert(cfg['default']['tests']))
    args.threads_no = args.threads_no or convert(cfg['default']['threads_no'], int)
    args.exec = args.exec or convert(cfg['default']['exec'], str)
    args.progress = args.progress or convert(cfg['default']['progress'], bool)
    args.benchmark = args.benchmark or convert(cfg['default']['benchmark'])
    args.tester_arguments = args.tester_arguments or cfg['default']['tester_arguments'] 
    args.data = args.data or cfg['default']['data']
    args.scale = args.scale or convert(cfg['default']['scale'], float)
    args.scoring = args.scoring or convert(cfg['default']['scoring'])
    args.sorting = args.sorting or convert(cfg['default']['sorting'])
    
    # Mode: Find
    if args.find:
        assert len(args.find) in [2, 3]
        assert args.find[1][-1] in ['-', '+']
        results = load_res_file(args.find[0] + cfg['general']['results_ext'])
        tests = results.keys()
        for filter in args.filters or []:
            tests = apply_filter(tests, results, filter)
        var = args.find[1][:-1]
        ascending = args.find[1][-1] == '+'
        
        ordered_tests = [test for _, test in sorted(zip([results[test][var] for test in tests], tests), reverse=not ascending)]
        
        if len(args.find) == 3:
            ordered_tests = ordered_tests[:int(args.find[2])]
        
        for test in ordered_tests:
            print(json.dumps(results[test]))
        sys.exit(0)
        
    # Parse args.tests
    if args.tests is None:
        pass
    elif isinstance(args.tests, int):
        args.tests = list(range(1, args.tests + 1))
    elif re.search('[a-zA-Z]', args.tests):
        if not os.path.exists(args.tests):
            print(f'Fatal Error: Cannot locate {args.tests} file')
            sys.exit(1)

        with open(args.tests) as f:
            lines = f.read().splitlines()
        assert len(lines) > 0
        
        if isinstance(try_str_to_numeric(lines[0]), int):
            args.tests = [int(line) for line in lines]
        else:
            args.tests = [json.loads(line)['id'] for line in lines]
    else:
        assert '-' in args.tests
        lo, hi = args.tests.split('-')
        lo = try_str_to_numeric(lo)
        hi = try_str_to_numeric(hi)
        args.tests = list(range(lo, hi + 1))
    
    # Mode: Summary 
    if args.show:
        results_files = find_res_files(cfg['general']['results_dir'])
        if args.sorting == 'name':
            results_files = sorted(results_files)
        elif args.sorting == 'date':
            results_files = [result_file for _, result_file in sorted(zip([os.path.getmtime(result_file) for result_file in results_files], results_files))]
            
        results = {os.path.basename(file).split('.')[0]: load_res_file(file) for file in results_files}
        
        if args.data == 'LATEST':
            _, args.data = sorted(zip([os.path.getmtime(result_file) for result_file in results_files], results_files))[-1]
        data_file = load_res_file(args.data) if args.data and os.path.isfile(args.data) else None
        
        show_summary(results, tests=args.tests, data=data_file, groups=args.groups, filters=args.filters)
        sys.exit(0)


    # Mode: Run tests
    if not os.path.exists(cfg['general']['tests_dir']):
        os.mkdir(cfg['general']['tests_dir'])
        
    if not args.tests:
        print('[Fatal Error] You need to specify tests to run, use --tests option')
        sys.exit(1)
    
    assert args.threads_no >= 1
    fout = sys.stdout
    if args.name:
        fout = open(f'{cfg["general"]["results_dir"]}/{args.name}{cfg["general"]["results_ext"]}', 'w')
        
    #TODO: add errors handling/warnings for benchmark file (file not existing, no full test coverage)
    benchmark = load_res_file(args.benchmark + cfg['general']['results_ext']) if args.benchmark else None
    
    try:
        start_time = time.time()
        for id in args.tests:
            tests_queue.put(id)
        tests_left = args.tests
        workers = [Thread(target=worker) for _ in range(args.threads_no)]
        for worker in workers:
            worker.start()
        
        sum_scores = 0
        log_scores = 0
        benchmark_log_scores = 0
        results = {}
        processed = 0
        
        while tests_left:
            result = results_queue.get()
            results[result['id']] = result
            assert result['id'] in tests_left
            while tests_left and tests_left[0] in results:
                processed += 1
                seed = tests_left[0]
                print(json.dumps(results[seed]), file=fout, flush=True)
                tests_left = tests_left[1:]
                sum_scores += results[seed]['score'] if results[seed]['score'] > 0 else 0
                log_scores += math.log(results[seed]['score']) if results[seed]['score'] > 0 else 0
                if args.progress and args.name:
                    output = f'Progress: {processed} / {processed+len(tests_left)}   Time: {time.time() - start_time : .3f}'
                    if args.benchmark:
                        benchmark_log_scores += math.log(benchmark[seed]['score'] if benchmark[seed]['score'] > 0 else 0)
                        output += f'   Scores: {log_scores / processed : .6f} vs {benchmark_log_scores / processed : .6f}'
                    print(f'\r{output}                       ', end='', file=sys.stderr)
                    sys.stderr.flush()
                    time.sleep(0.001)
    except KeyboardInterrupt:
        print('\nInterrupted by user', file=sys.stderr)
        os._exit(1)
        
    print(file=sys.stderr)
        
    print("Time:", time.time() - start_time, file=sys.stderr)
    print("Avg Score:", sum_scores / len(results), file=sys.stderr)
    print("Avg Log Scores:", log_scores / len(results), file=sys.stderr)
    
    