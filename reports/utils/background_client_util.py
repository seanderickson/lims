##############################################################################
# Command line background job processing utility:
# - services a background job from the command line using the 
# reports.utils.background_processor.BackgroundClient
# - the BackgroundClient will utilize the Django server in "headless" mode
# to perform API requests.
##############################################################################

import argparse
import datetime
import logging
import multiprocessing
import os.path
import socket
import subprocess

import django
from django.conf import settings

from reports import InformationError, ValidationError
from reports.utils import parse_credentials
import reports.utils.background_processor


logger = logging.getLogger(__name__)


def execute_from_python(job_id, sbatch=False):
    '''
    Utility method to invoke from the running server.
    
    @see settings.BACKGROUND_PROCESSOR
    
    @param sbatch if true, requires "sbatch_settings" in the 
    BACKGROUND_PROCESSOR settings
    @param keep_stdout (for testing) set True to use STDOUT (for non-sbatch only)
    '''
    logger.info('using settings.BACKGROUND_PROCESSOR: %r', 
        settings.BACKGROUND_PROCESSOR)
    
    check_settings = set(['post_data_directory','job_output_directory',
        'credential_file', 'python_environ_script', 'background_process_script'])
    
    if not check_settings.issubset(
        set(settings.BACKGROUND_PROCESSOR.keys())):
        raise ValidationError(
            key='settings.BACKGROUND_PROCESSOR', 
            msg='missing required entries: %s' 
                % (check_settings-set(settings.BACKGROUND_PROCESSOR.keys())))
    
    job_output_dir = settings.BACKGROUND_PROCESSOR['job_output_directory']
    if not os.path.exists(job_output_dir):
        os.makedirs(job_output_dir)
    credential_file = settings.BACKGROUND_PROCESSOR['credential_file']
    python_environ_script = settings.BACKGROUND_PROCESSOR['python_environ_script']

    # background_process_script = settings.BACKGROUND_PROCESSOR['background_process_script']
    background_process_script = os.path.abspath(__file__)
    logger.info('this file: %r', background_process_script)
    
    output_stdout = '%d.stdout'%job_id
    output_stdout = os.path.abspath(os.path.join(job_output_dir,output_stdout))
    output_stderr = '%d.stderr'%job_id
    output_stderr = os.path.abspath(os.path.join(job_output_dir,output_stderr))

    run_sh_args = [
        python_environ_script, background_process_script, 
        '--job_id', str(job_id), '--c', credential_file, '-vv']
    full_args = []
    
    if sbatch is True:
        os.putenv('USER', 'sde4')
        full_args.append('/usr/local/bin/sbatch')

        sbatch_settings = settings.BACKGROUND_PROCESSOR.get('sbatch_settings')
        if sbatch_settings is None:
            raise InformationError(
                key='sbatch_settings', 
                msg='missing from the BACKGROUND_PROCESSOR settings')
        
        sbatch_settings['output'] = output_stdout
        sbatch_settings['error'] = output_stderr
        sbatch_settings['job-name'] = 'ss_{}'.format(job_id)
        sbatch_args = []
        for k,v in sbatch_settings.items():
            sbatch_args.extend(['--%s=%s' % (k, '%s'%str(v))])
        full_args.extend(sbatch_args)
        full_args.append('-vvv')

    full_args.extend(run_sh_args)
    
    logger.info('full args: %r', full_args)
    
    if sbatch is True:
        logger.info('sbatch specified, invoke sbatch and wait for output...')
        logger.info('full command %s: ', ' '.join(full_args))
        try:
            output = \
                subprocess.check_output(full_args, stderr=subprocess.STDOUT)
            logger.info('ran, output: %r', output)
            # TODO: parse the SLURM process ID from the output
            return output
        except subprocess.CalledProcessError, e:
            logger.error('subprocess.CalledProcessError: output: %r', e.output)
            raise
    else:
        # NOTE: Several shell options are presented:
        # - shell operation is a simple drop in replacement for sbatch
        # NOTE: running as a multiprocess process should simplify this, but
        # properly routing the logs for multiprocessing requires rework
        
        hostname = socket.gethostname()
        logger.info('sbatch not specified, run directly in the shell on %r', hostname)
        
        SHELL_ASYNC_EXEC = 'SHELL_ASYNC_EXEC'
        SHELL_ASYNC_SUBSHELL = 'SHELL_ASYNC_SUBSHELL'
        SHELL_SYNC_PIPES = 'SHELL_SYNCHRONOUS_TO_PIPES'
        
        type = SHELL_ASYNC_EXEC
        logger.info('run as %r...', type)
        
        # NOTE: all options spawn child processes, so the PID is the parent
        if type == SHELL_ASYNC_EXEC:
            logger.info('stdout: %r', output_stdout)
            stdout = open(output_stdout, "w")
            logger.info('stderr: %r', output_stderr)
            stderr = open(output_stderr, "w")

            logger.info('full_args: %r', full_args)
            process = subprocess.Popen(full_args, stdout=stdout, stderr=stderr)
            logger.info('ran as %r, process id: %r', type, process.pid)
            return process.pid
        elif type == SHELL_ASYNC_SUBSHELL:
            full_args.append('>%s'%output_stdout)
            full_args.append('2>%s'%output_stderr)
            full_args = ' '.join(full_args)    
            process = subprocess.Popen(full_args, shell=True)
            logger.info('ran as %r, process id: %r', type, process.pid)
            return process.pid
        elif type == SHELL_SYNC_PIPES:
            # NOTE: testing only, job will block before returning
            logger.info('full_args: %r', full_args)
            process = subprocess.Popen(full_args, 
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            logger.info('ran async process id: %r', process.pid)

            outs, errs = process.communicate() 
            print 'STDOUT', outs
            print 'STDERR', errs
            return process.pid
        else:
            raise Exception('unknown run type: %r', type)

parser = argparse.ArgumentParser(description=
    'background_client_util: '
    'Service a pending background job from the command line using the ' 
    'reports.utils.background_processor.BackgroundClient '
    'to access the Django server in "headless" mode to perform API requests. '
    'NOTE: requires a valid DJANGO_SETTINGS_MODULE in the environment.')

# Admin credentials
parser.add_argument(
    '-U', '--username',
    help='username for the api authentication')
parser.add_argument(
    '-p', '--password',
    help='user password for the api authentication')
parser.add_argument(
    '-c', '--credential_file',
    help='credential file containing the username:password for api '
    'authentication')

parser.add_argument(
    '-j', '--job_id', required=True, type=int, 
    help='Job ID to process')

parser.add_argument(
    '-run_sbatch', '--run_sbatch', action='store_true',
    help='calls this script using sbatch settings (for testing)'
        'Note: either "run_sbatch" or "run_subprocess" or "run_multiprocess" '
        'may be specified')
parser.add_argument(
    '-run_subprocess', '--run_subprocess', action='store_true',
    help='calls this script using a subprocess (for testing). '
        'Note: either "run_sbatch" or "run_subprocess" or "run_multiprocess" '
        'may be specified')
parser.add_argument(
    '-run_multiprocess', '--run_multiprocess', action='store_true',
    help='calls this script using a multiprocessing.process (for testing). '
        'Note: either "run_sbatch" or "run_subprocess" or "run_multiprocess" '
        'may be specified')


parser.add_argument(
    '-v', '--verbose', dest='verbose', action='count',
    help="Increase verbosity (specify multiple times for more)")    
    
        
if __name__ == "__main__":
    args = parser.parse_args()
    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
        DEBUG=True
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        

    if args.credential_file:
        username,password = parse_credentials(args.credential_file)
    if username is None:
        username = args.username
        if username is None:
            parser.error(
                'username is required if not specifying the credential_file')
        password = args.password
        if not password:
            password = getpass.getpass()

    if args.run_sbatch is True:
        print 'Send the job to sbatch...', args.job_id
        output = execute_from_python(args.job_id, sbatch=True)
        print 'job output:', output
    elif args.run_subprocess is True:
        print 'Processing the job as a subprocess call...', args.job_id
        output = execute_from_python(args.job_id, sbatch=False)
        print 'job output:', output
    elif args.run_multiprocess is True:
        print 'Processing the job as a multiprocessing.Process...', args.job_id
        django.setup()
        api_client = \
            reports.utils.background_processor.ApiClient(username, password)
        background_client = \
            reports.utils.background_processor.BackgroundClient(api_client)
        p = multiprocessing.Process(
            target=background_client.service,args=(args.job_id,) )
        p.start();
        logger.info('started multiprocessing.Process: %r', p.pid)
    else:
        django.setup()
        api_client = \
            reports.utils.background_processor.ApiClient(username, password)
        background_client = \
            reports.utils.background_processor.BackgroundClient(api_client)
        response = background_client.service(args.job_id)    

    print 'exit background processing service', datetime.datetime.now()
    