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
    #general.add_argument('-refEnv', '--refEnviro', type = str, help = 'Path to species per Environment reference')
    general.add_argument( '-ref','--ref', type = str, help = 'Path to SILVA alignment reference')
    general.add_argument( '-refTax','--refTax', type = str, help = 'Path to SILVA alignment TAX reference')
    args = general.parse_args()
    return(args)



def main(args):
    make_install_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    DB = '/'.join(make_install_home.split('/')[:-1] + ['DB'])
    #try:
    #    os.mkdir(DB)
    #except OSError as e: #[Errno 17] File exists: 'outputDir'
    #    if e.errno != 17:
    #        raise
    #    elif not args.force_overwrite:
    #        print('warning: The directory {} already exists. Please, remove it or choose other name for the output directory'.format(DB))
    #        exit(-17)
    #    else:
    #        shutil.rmtree(DB, ignore_errors = True)
    #        os.mkdir(DB)
    os.chdir(DB)
    
    if (not args.refTax) or (not args.ref):
        silva_version = 'silva.nr_v138.tgz'
        silva_url = 'https://mothur.s3.us-east-2.amazonaws.com/wiki/' + silva_version
        print('Please wait until the downloading has completely finished')
        silva = requests.get(silva_url, stream = True)
        #silva.status_code == requests.codes.ok # Be sure that the request is ok!
        silva.raise_for_status() # check the url is ok!
        with open(silva_version, 'wb') as f: # write the content to a file
            f.write(silva.raw.read())
        with tarfile.open(silva_version, 'r') as tar_ref:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar_ref)
    else:
       os.symlink(args.ref, '{}/{}'.format(DB, args.ref.split('/')[-1]))
       os.symlink(args.refTax, '{}/{}'.format(DB, args.refTax.split('/')[-1]))
 
        
    #if not args.refEnviro:
        # Species per Enviro: always download!
    #    env_url = 'https://github.com/ggnatalia/MMs/tree/main/DB' # Change for raw link
    #    env = requests.get(env_url)
    #    env.raise_for_status()
    #else:
    #    os.symlink(args.refEnviro, '.')   
    
    print('Your DB are in:\n{}\n'.format(DB))


################################################################################################################
if __name__ == '__main__':
    main(parse_arguments())
