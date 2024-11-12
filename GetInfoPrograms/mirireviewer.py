from urllib.request import urlopen
from bs4 import BeautifulSoup

def get_reviewer(program):
    print(f'Searching for PID {program} info...')
    url = f"https://www.stsci.edu/cgi-bin/get-proposal-info?id={program}&observatory=JWST"
    html = urlopen(url).read()
    soup = BeautifulSoup(html, features="html.parser")
    paragraphs = soup.find_all('p')
    if 'MIRI Reviewer' in soup.get_text():
        for p in paragraphs:
            line = p.get_text()
            if 'MIRI Reviewer' in line:
                if '@' in " ".join(line.split(" ")):
                    idx = [i for i in range(len(line.split(" "))) if '@' in line.split(" ")[i]][0]
                    reviewer = " ".join(line.split(" ")[2:idx])
                else:
                    reviewer = " ".join(line.split(" ")[2:-1])
    else:
        reviewer = "N/A"
    return reviewer