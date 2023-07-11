import requests, sys, json, argparse, re, os.path, os

BASE_URL = 'api.icfpcontest.com'

HEADER = {'authorization' : 'Bearer TOKEN'}

def extract_results(path):
    path = os.path.basename(path)
    match = re.match(r'test(?P<id>\d+)_(?P<score>\d+)', path)
    if match is None:
        return None, None
    return int(match['id']), int(match['score'])


def submit(path):
    print('Submitting:', path)
    id, score = extract_results(path)
    assert id is not None and score is not None
    data = open(path, 'r').read()
    j = requests.post('https://api.icfpcontest.com/submission', json={"contents": data, "problem_id": id}, headers=HEADER).json()
    if 'Failure' in j: 
        print(j)
        sys.exit(1)


def get_scores():
    print('Downloading Userboard')
    j = requests.get(f'https://{BASE_URL}/userboard', headers=HEADER).json()
    if 'Failure' in j: 
        print(j)
        sys.exit(1)
    return j['Success']['problems']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ICFP 2023 API; for help on a specific mode, type cli.py <mode> -h')
    parser.set_defaults(mode='unknown')
    subparsers = parser.add_subparsers(title='modes')

    parser_download = subparsers.add_parser('download', aliases=['d'], help='download problem')
    parser_download.set_defaults(mode='download')
    parser_download.add_argument('id', type=str, help='number or range')
    parser_download.add_argument('dir', default='.', type=str, nargs='?', help='directory')

    parser_submit = subparsers.add_parser('submit', aliases=['s'], help='submit solution')
    parser_submit.set_defaults(mode='submit')
    parser_submit.add_argument('path', type=str, help='test to submit')

    parser_submit_all = subparsers.add_parser('submit_all', aliases=['sa'], help='submit all solutions that improve best scores')
    parser_submit_all.set_defaults(mode='submit_all')
    parser_submit_all.add_argument('dir', default='.', type=str, nargs='?', help='directory')

    parser_userboard = subparsers.add_parser('userboard', aliases=['u'], help='get userboard')
    parser_userboard.set_defaults(mode='userboard')

    args = parser.parse_args()
    

    if args.mode == 'unknown':
        parser.print_help()

    elif args.mode == 'download':
        ids = [int(args.id)] if '-' not in args.id else list(range(int(args.id.split('-')[0]), int(args.id.split('-')[1])+1))
        for id in ids:
            print('Downloading', id)
            test = requests.get(f'http://cdn.icfpcontest.com/problems/{id}.json').json()
            with open(f'{args.dir}/test{id}.json', 'w') as f:
                json.dump(test, f)

    elif args.mode == 'submit':
        submit(args.path)

    elif args.mode == 'submit_all':
        scores = get_scores()
        scores = [int(score) if score is not None else 0 for score in scores]
        best_results = {}
        for path in os.listdir(args.dir):
            id, score = extract_results(path)
            if id is None or score is None:
                continue
            if id not in best_results or score > best_results[id][0]:
                best_results[id] = (score, path)

        total_improvement = 0
        submits = []
        for id, (score, path) in best_results.items():
            if score > scores[id-1]:
                print(f'Test #{id} improves score from {scores[id-1]} to {score}: +{score - scores[id-1]}')
                total_improvement += score - scores[id-1]
                submits.append(path)

        print('Total Improvement:', total_improvement)
        for path in submits:
            submit(path)

    elif args.mode == 'userboard':
        scores = get_scores()
        n_tests = len(scores)
        sum = 0
        sumP1 = 0
        sumP2 = 0

        for id, score in enumerate(scores):
            print(f'#{id+1}\t: {int(score) if score is not None else "N/A":>12}')
            score = int(score or 0)
            sum += score
            if id+1 <= 55:
                sumP1 += score
            else:
                sumP2 += score

        print('P1:', sumP1)
        print('P2:', sumP2)
        print('Total:', sum)
