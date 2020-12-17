#!/usr/bin/env python3

import argparse
import os
import shutil
import requests
import tarfile
# Script for downloading SILVA db and the species per Environment table.


def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    general.add_argument('-refEnv', '--refEnviro', type = str, help = 'Path to species per Environment reference')
    general.add_argument( '-ref','--ref', type = str, help = 'Path to SILVA alignment reference')
    general.add_argument( '-refTax','--refTax', type = str, help = 'Path to SILVA alignment TAX reference')
    general.add_argument('--force-overwrite', action = 'store_true', help = 'Force overwrite if the output directory already exists')
    args = general.parse_args()
    return(args)



def main(args):
    make_install_home = abspath(dirname(realpath(__file__)))
    DB = '/'.join(abspath(dirname(realpath(__file__))).split('/')[:-1] + ['DB'])
    try:
        os.mkdir(DB)
    except OSError as e: #[Errno 17] File exists: 'outputDir'
        if e.errno != 17:
            raise
        elif not args.force_overwrite:
            print('warning: The directory {} already exists. Please, remove it or choose other name for the output directory'.format(DB))
            exit(-17)
        else:
            shutil.rmtree(DB, ignore_errors = True)
            os.mkdir(DB)
    os.chdir(DB)
    
    if (not args.refTax) and (not args.ref):
        silva_version = 'silva.nr_v138.tgz'
        silva_url = 'https://mothur.s3.us-east-2.amazonaws.com/wiki/' + silva_version
        silva = requests.get(silva_url, stream = True)
        #silva.status_code == requests.codes.ok # Be sure that the request is ok!
        silva.raise_for_status() # check the url is ok!
        with open(silva_version, 'wb') as f: # write the content to a file
            f.write(silva.raw.read())
        with tarfile.open(silva_version, 'r') as tar_ref:
            tar_ref.extractall()
    else:
        os.symlink(args.ref, '.')
        os.symlink(args.refTax, '.')
        
        
    if not args.refEnviro:
        # Species per Enviro: always download!
        env_url = 'https://raw.githubusercontent.com/ggnatalia/MMs/main/DB/SpeciesperEnviro.tsv?token=AMC24ROKCNQHCPXIMI7BWTC73NXVK'
        env = requests.get(env_url)
        env.raise_for_status()
    else:
        os.symlink(args.refEnviro, '.')   
    
    print('Your DB are in:\n{}\n'.format(DB))


################################################################################################################
if __name__ == '__main__':
    main(parse_arguments())
