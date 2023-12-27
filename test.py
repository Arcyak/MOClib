import astropy_healpix

import subprocess
import sys
import re


def extract_votable(catalogue: str, out_max: int = 10000):
    """get VOTable from a VizieR catalogue for further studies.
    :param catalogue : VizieR catalogue name
    :param out_max: max number of objects, default: 10000
    """

    result = []
    proc = subprocess.Popen(["/bin/sh"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)

    url = f'https://vizier.cds.unistra.fr/viz-bin/votable?-source={catalogue}&-out.max={out_max}'

    # execute the command
    command = f"curl {url}"
    proc.stdin.write(bytes(command, "ascii"))  # send command in STDIN
    proc.stdin.close()  # close STDIN

    # read data (STDOUT)
    for line in proc.stdout:
        line = line.decode().lstrip()
        record = [value.strip() for value in re.split(r" +", line)]
        result.append(record)
    proc.stdout.close()
    print(result)

    # read if errors
    for line in proc.stderr:
        sys.stderr.write(f"(error) {line.decode()}\n")
    proc.stderr.close()

    return result


if __name__ == '__main__':
    print(extract_votable('II/7A/catalog'))

