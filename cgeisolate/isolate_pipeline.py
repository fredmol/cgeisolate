import os
import sys
import subprocess

from cgeisolate import kma

def isolate_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    kma.KMARunner(args.input,
              args.output + "/reference_mapping",
              args.db_dir + '/bac_db',
              "-ID 50 -nf -mem_mode -sasm -ef -1t1").run()
    return 'isolate_pipeline'