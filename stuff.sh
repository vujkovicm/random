# print duplicate lines
awk 'seen[$0]++' filename 
