from pprint import pprint
import requests

URL = 'https://discosweb.esoc.esa.int'
token = 'ImI5MjBiMzM2LTk2ZjktNDA0ZC1hZGFiLTNjMmFiZjc4NGFkYSI.0_k20tZIVHaaq-89tngta3rmZd0'

response = requests.get(
    f'{URL}/api/launches',
    headers={
        'Authorization': f'Bearer {token}',
        'DiscosWeb-Api-Version': '2',
    },
    params={
        'filter': "ge(epoch,epoch:'2010-01-01')&le(epoch,epoch:'2017-12-31')",
        'sort': 'epoch',
        'page[size]': '100',
        'page[number]': '1'
    },
)
doc = response.json()
if response.ok:
    pprint(doc['data'])
else:
    pprint(doc['errors'])
print(response.url)