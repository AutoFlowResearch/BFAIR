from isatools.io import isatab_parser as ip

"""investigation = Investigation()
investigation.studies.append(Study())  # adds a new default Study object to investigation


print(investigation.studies)"""

rec = ip.parse("../sample_data")

print(rec.studies[0])


