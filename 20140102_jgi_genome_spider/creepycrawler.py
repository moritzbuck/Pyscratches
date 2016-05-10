""" small script downloading genomes from jgi

you need an account on jgi and a "project list"-file which  can be downloaded from jgi by using the export to excel button (with the inclure ressources box checked) on the project list page (http://genome.jgi.doe.gov/genome-projects/)

Usage:
  creepycrawler.py [-o <dir>  -p <password> -u <username> -l <MB>] (-i <file>)

Options:
   -h --help       Show this screen.
   -o <dir>        target
   -i <file>       project list
   -p <password>   jgi password
   -u <username>   jgi username
   -l <MB>         file size limit in MB

"""
from docopt import docopt
import pycurl
import StringIO
from pandas.io import parsers
import elementtree.ElementTree as ET
import re
import os
import sys

def mkfolder(element,upfolder):
    upfolder=upfolder+element.attrib['name']+'/'
    if not os.path.exists(upfolder):
        os.mkdir(upfolder)
    for folder in element.findall('folder'):
        mkfolder(folder,upfolder)
    for fil in element.findall('file'):
        mkfile(fil,upfolder)

def mkfile(element, upfolder):
    url = element.attrib['url']
    chunks = url.split('/')
    target = 'http://genome.jgi.doe.gov' + url
    path = upfolder + chunks[1] + chunks[-1]
    
    if not os.path.exists(path):
        c.setopt(pycurl.URL, target)
        
        with open(path, 'wb') as out:
            c.setopt(pycurl.WRITEFUNCTION, out.write)
            c.perform()


if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    user = arguments['-u']
    psswd = arguments['-p']
else :
    user = "mrtz.buck@gmail.com"
    psswd = "9Y4xQ2dh."

if arguments['-l'] == None:
    limit=1024
else:
    limit=int(arguments['-l'])

print(limit)

#sys.exit()
    
c = pycurl.Curl()
c.setopt(pycurl.URL,"https://signon.jgi.doe.gov/signon/create")
c.setopt(c.COOKIEFILE, '')
c.setopt(pycurl.POST, 1)
c.setopt(pycurl.HTTPPOST, [('login', user), ('password', psswd)])
c.perform()

c.setopt(pycurl.POST, 0) 


url = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=KeggMap&page=keggMapRelated&map_id=map02010&taxon_oid=2626542262"

c.setopt(pycurl.HEADER,1)

c.setopt(pycurl.URL,url)
out = StringIO.StringIO()
c.setopt(pycurl.WRITEFUNCTION, out.write)
c.perform()

temp = [l for l in out.getvalue().split("\n") if "<img" in l and "map" in l ][0]
map = re.search("src='(.+?)'", temp).groups()[0]
