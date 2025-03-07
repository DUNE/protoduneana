from argparse import ArgumentParser as ap
import subprocess
import yaml
import os.path

def check_fcls(fits):
  bad_files = []
  for fit, outdir, mc_file, data_file, N in fits:
    if not os.path.isfile(fit):
      print('Found bad fcl file', fit)
      bad_files.append(fit)

  if len(bad_files) > 0:
    raise ValueError('Bad fcl files', bad_files)
  else:
    print('All fcls good')

def check_input_files(fits, mc_dir, data_dir):
  bad_mc_files = []
  bad_data_files = []

  all_mc_files = []
  all_data_files = []
  for fit, outdir, mc_file, data_file, N in fits:

    #full_mc_file = (mc_file if '/pnfs' in mc_file else mc_dir + '/' + mc_file)
    #full_data_file = (data_file if '/pnfs' in data_file else data_dir + '/' + data_file)
    full_mc_file = make_file(mc_file, mc_dir)
    full_data_file = make_file(data_file, data_dir)

    all_mc_files.append(full_mc_file)
    all_data_files.append(full_data_file)

  all_mc_files = set(all_mc_files)
  all_data_files = set(all_data_files)
  for f in all_mc_files:
    if not os.path.isfile(f):
      bad_mc_files.append(f)
  for f in all_data_files:
    if not os.path.isfile(f):
      bad_data_files.append(f)

  if len(bad_mc_files) > 0 or len(bad_data_files) > 0:
    raise ValueError('Bad mc files:', bad_mc_files,
                     '\nBad data files:', bad_data_files)
  else:
    print('All input files good')

def make_file(f, fdir):
  return (f if '/pnfs' in f else fdir + '/' + f)
  

if __name__ == '__main__':

  parser = ap()
  parser.add_argument('-i', type=str)
  parser.add_argument('--mid_dir', type=str, default='')
  parser.add_argument('--memory', type=str, default='13999MB')
  parser.add_argument('--lifetime', type=str, default='3h')
  parser.add_argument('--blacklist', type=str, default=['RAL'], nargs='+')
  parser.add_argument('--pduneana_tar', type=str, default=None)
  parser.add_argument('--dry_run', action='store_true')
  parser.add_argument('--extra_dry_run', action='store_true')
  parser.add_argument('--config', type=str, default=None)
  args = parser.parse_args()

  with open(args.i, 'r') as f: config = yaml.safe_load(f)

  topdir = config['topdir']

  fits = config['fits']
  check_fcls(fits)

  mc_dir = config['mc_dir']
  data_dir = config['data_dir']
  check_input_files(fits, mc_dir, data_dir)

  for fit, outdir, mc_file, data_file, N in fits:
    print(fit, outdir, N)

    cmd = [
      'python',
      'submit_fits_local.py', '--type', fit.replace('.fcl', ''),
      '--input_file', make_file(mc_file, mc_dir),
      '--data_input', make_file(data_file, data_dir),
      '--memory', args.memory,
      '--lifetime', args.lifetime,
    ]

    if args.pduneana_tar is not None:
      if not os.path.isfile(args.pduneana_tar):
        raise Exception(args.pduneana_tar, 'not found')
      cmd += ['--pduneana_tar', args.pduneana_tar]

    if args.blacklist is not None and len(args.blacklist) > 0:
      cmd += ['--blacklist'] + [b for b in args.blacklist]
    
    cmd += [
      '--output_dir', f'{topdir}/{args.mid_dir}/{outdir}',
      '--alloutput',
    ]

    if args.config is not None:
      cmd += ['--config', args.config]

    if N > 1:
      cmd += ['--multiple', '-N', str(N)]

    if args.dry_run: cmd += ['--dry_run']

    if args.extra_dry_run:
      print('Will run')
      print(cmd[0], '\\')
      for i in range(1, len(cmd), 2):
        if i < len(cmd)-1:
          print('\t',cmd[i], cmd[i+1], '\\')
        else:
          print('\t',cmd[i])
    else:
      subprocess.run(cmd) 
