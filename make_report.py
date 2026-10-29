"""Build Project_Report.docx from README.md + SAP.md (needs pandoc). python3 make_report.py"""
import re, subprocess

readme = open("README.md").read()
readme = re.sub(r"^!\[tests\].*\n\n?", "", readme, flags=re.M)              # remote badge: not for a Word doc
readme = re.sub(r"^\*\*Stack:\*\*.*\n.*\n", "", readme, flags=re.M)          # link line points at this file
readme = readme.split("## Repo map")[0]
readme = readme.replace("[statistical analysis plan](SAP.md)", "statistical analysis plan (Appendix)")
readme = readme.replace("([SAP.md](SAP.md))", "(Appendix)")
sap = open("SAP.md").read().replace("# Statistical analysis plan", "# Appendix: statistical analysis plan")
doc = ("---\ntitle: Survival extrapolation and network meta-analysis for HTA\nsubtitle: Project report, Sharon Okolo\n"
       "---\n\nCode, data outputs and the Excel model: https://github.com/SharonOkolo/hta-survival-nma\n\n"
       + readme.split("\n", 1)[1] + "\n\n" + sap)
open("/tmp/hta_report.md", "w").write(doc)
subprocess.run(["pandoc", "/tmp/hta_report.md", "-o", "Project_Report.docx", "--resource-path=.",
                "--reference-doc=../../templates/cv_reference.docx"], check=True)
print("wrote Project_Report.docx")
