import requests
import tarfile
import time
import os
import shutil
import random
from tqdm import tqdm
import logging
from Bio import SeqIO
logger = logging.getLogger(__name__)

# Adapted from https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

host_url="https://api.colabfold.com"
user_agent = "ProteinCartography/0.4 (Arcadia Science) python-requests/2.0.1"
headers = {}
headers['User-Agent'] = user_agent
submission_endpoint = "ticket/msa"

def run_mmseqs2(query, output_folder, output_file, protid):

  def submit(seqs, mode, N=101):
    n, query = N, ""
    for seq in seqs:
      query += f">{n}\n{seq}\n"
      n += 1

    while True:
      error_count = 0
      try:
        # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
        # "good practice to set connect timeouts to slightly larger than a multiple of 3"
        res = requests.post(f'{host_url}/{submission_endpoint}', data={ 'q': query, 'mode': mode }, timeout=6.02, headers=headers)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while submitting to MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break

    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def status(ID):
    while True:
      error_count = 0
      try:
        res = requests.get(f'{host_url}/ticket/{ID}', timeout=6.02, headers=headers)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while fetching status from MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def download(ID, path):
    error_count = 0
    while True:
      try:
        res = requests.get(f'{host_url}/result/download/{ID}', timeout=6.02, headers=headers)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while fetching result from MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break
    with open(path,"wb") as out: out.write(res.content)

  # process input x
  seqs = [str(SeqIO.read(query, "fasta").seq)]

  # setup mode
  mode = "nofilter"

  # define path
  if not os.path.isdir(output_folder): os.mkdir(output_folder)

  # call mmseqs2 api
  tar_gz_file = f'{output_folder}/out.tar.gz'
  N,REDO = 101,True

  # deduplicate and keep track of order
  seqs_unique = []
  #TODO this might be slow for large sets
  [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
  Ms = [N + seqs_unique.index(seq) for seq in seqs]
  # lets do it!
  if not os.path.isfile(tar_gz_file):
    TIME_ESTIMATE = 150 * len(seqs_unique)
    with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
      while REDO:
        pbar.set_description("SUBMIT")

        # Resubmit job until it goes through
        out = submit(seqs_unique, mode, N)
        while out["status"] in ["UNKNOWN", "RATELIMIT"]:
          sleep_time = 5 + random.randint(0, 5)
          logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
          # resubmit
          time.sleep(sleep_time)
          out = submit(seqs_unique, mode, N)

        if out["status"] == "ERROR":
          raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

        if out["status"] == "MAINTENANCE":
          raise Exception(f'MMseqs2 API is undergoing maintenance. Please try again in a few minutes.')

        # wait for job to finish
        ID,TIME = out["id"],0
        pbar.set_description(out["status"])
        while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
          t = 5 + random.randint(0,5)
          logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
          time.sleep(t)
          out = status(ID)
          pbar.set_description(out["status"])
          if out["status"] == "RUNNING":
            TIME += t
            pbar.update(n=t)
          #if TIME > 900 and out["status"] != "COMPLETE":
          #  # something failed on the server side, need to resubmit
          #  N += 1
          #  break

        if out["status"] == "COMPLETE":
          if TIME < TIME_ESTIMATE:
            pbar.update(n=(TIME_ESTIMATE-TIME))
          REDO = False

        if out["status"] == "ERROR":
          REDO = False
          raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

      # Download results
      download(ID, tar_gz_file)


  # extract a3m files
  with tarfile.open(tar_gz_file) as tar_gz:
    for member in tar_gz.getmembers():
      if member.name == "uniref.a3m":
        tar_gz.extract(member, f"{output_folder}/{protid}")

  os.rename(f"{output_folder}/{protid}/uniref.a3m", f"{output_folder}/{protid}/{protid}_mmseqs2_hits.a3m")
  shutil.move(f"{output_folder}/{protid}/{protid}_mmseqs2_hits.a3m", output_file)
  os.rmdir(f"{output_folder}/{protid}")
  os.remove(tar_gz_file)