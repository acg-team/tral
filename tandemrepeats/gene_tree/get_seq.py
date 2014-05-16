import httplib2, sys
 
http = httplib2.Http(".cache")
  
server = "http://test.rest.ensemblgenomes.org"
ext = "/sequence/id/ENSPVAG00000014049"
resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"text/fasta"})
   
if not resp.status == 200:
    print("Invalid response: ", resp.status)
    sys.exit()       

print(content)
