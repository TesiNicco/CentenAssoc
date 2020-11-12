#!/usr/bin/python
# PACKAGES
import re
import werkzeug
werkzeug.cached_property = werkzeug.utils.cached_property
from robobrowser import RoboBrowser
import sys

# MAIN
## read inputs go and pvalues
d = open("RESULTS/revigo_inp.txt").read()

## set connections
br = RoboBrowser(parser="lxml")
br.open("http://revigo.irb.hr/")

## manage parameters
clus_sim = float(sys.argv[1])
if clus_sim == 0.4:
    clus_sim = '0.40'
elif clus_sim == 0.5:
    clus_sim = '0.50'
elif clus_sim == 0.7:
    clus_sim = '0.70'

sim_meas = sys.argv[2]
if sim_meas == "Lin":
    sim_meas = 'Lin'
if sim_meas == "SimRel":
    sim_meas = 'SIMREL'

## fill in my data
form = br.get_form()
form["goList"].value = d
form["cutoff"].value = clus_sim
form["goSizes"].value = '9606'
form["measure"].value = sim_meas

## submit request
br.submit_form(form)

## get link to download stuff
download_csv_link = br.find("a", href=re.compile("export.jsp"))
br.follow_link(download_csv_link)
csv_content = br.response.content.decode("utf-8")

## write output
fout = open("RESULTS/revigo_out.csv", "w")
fout.write(csv_content)
fout.close()



