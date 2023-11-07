import os
import sys
import subprocess

from cgeisolate import kma

def isolate_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    kma.KMARunner(args.input,
              args.output + "/bacteria_alignment",
              args.db_dir + '/bac_db',
              "-ID 75 -md 5 -ont -1t1").run()

    kma.KMARunner(args.input,
              args.output + "/amr",
              args.db_dir + '/resfinder_db',
              "-ont -md 5").run()

    kma.KMARunner(args.input,
                  args.output + "/virulence",
                  args.db_dir + '/virulence_db',
                  "-ont -md 5").run()

    kma.KMARunner(args.input,
                  args.output + "/plasmid",
                  args.db_dir + '/plasmid_db',
                  "-ont -md 5").run()


    return 'isolate_pipeline'