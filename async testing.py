import aiohttp
import asyncio
import time

async def get_pokemon(session, url):
    async with session.post(url, headers={'Authorization' : 'Key 218fc836-c568-4e88-b12a-eaef0f538431'}, json={'amount' : 1, 'outcome' : 'YES', 'contractId' : '8b2EC2Au504RG4OCZnqd'}) as resp:
        pokemon = await resp.json(content_type=None)
        return pokemon["betId"]


async def main(length):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for number in range(1, length+1):
            url = 'https://manifold.markets/api/v0/bet'
            tasks.append(asyncio.ensure_future(get_pokemon(session, url)))
        original_pokemon = await asyncio.gather(*tasks)
        for pokemon in original_pokemon:
            print(pokemon)

start_time = time.time()
asyncio.run(main(5))
print('--- %s seconds ---' % (time.time() - start_time))


async def get_pokemon(session, url):
    async with session.post(url, headers={'authorization' : 'Bearer eyJhbGciOiJSUzI1NiIsImtpZCI6IjE2ZGE4NmU4MWJkNTllMGE4Y2YzNTgwNTJiYjUzYjUzYjE4MzA3NzMiLCJ0eXAiOiJKV1QifQ.eyJuYW1lIjoiTmVvbiBOdWtlIiwicGljdHVyZSI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FMbTV3dTJ5cXFuZHFNZzBOU0tfNlZSN0hFZ2UzeFBuQUpTOExHeVUwd0Y0dXc9czk2LWMiLCJpc3MiOiJodHRwczovL3NlY3VyZXRva2VuLmdvb2dsZS5jb20vbWFudGljLW1hcmtldHMiLCJhdWQiOiJtYW50aWMtbWFya2V0cyIsImF1dGhfdGltZSI6MTY3MzAxOTE0NCwidXNlcl9pZCI6ImdmVVJpaUdGWWVVb1V6dXpUdWJrTmFXZFN2TTIiLCJzdWIiOiJnZlVSaWlHRlllVW9VenV6VHVia05hV2RTdk0yIiwiaWF0IjoxNjgxODIwMzQwLCJleHAiOjE2ODE4MjM5NDAsImVtYWlsIjoidW5yYWh1bGJlYXRhYmxlQGdtYWlsLmNvbSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJmaXJlYmFzZSI6eyJpZGVudGl0aWVzIjp7Imdvb2dsZS5jb20iOlsiMTA4NjY4ODIyMjc4NDA5MDI3NjE0Il0sImVtYWlsIjpbInVucmFodWxiZWF0YWJsZUBnbWFpbC5jb20iXX0sInNpZ25faW5fcHJvdmlkZXIiOiJnb29nbGUuY29tIn19.RwYjut-bHcuIsvwvXO1BT7rahdhyM3t-r_fLKvsfDl4r3fiFv3oldJtB6xP8T9GTzE4Fu1EZmtggrom668lz-72wzlmw2RjotYCYRuvl4MgaycbEDjudCZK_L5uOR8rqjeNTiClK1BvtzDYpc7WcWr9KzW6iZ9EO6cVATJnWFO_UKldlMd9i81gRMR8pnt3t7MMXC-EmVdGZArybYVZ-G7kN4CGmXvsOzgcn66ZqGRpDwooxq8Fz8TY1pyj16j8lEAG32MC74Km_5QRhz550hK4OxzclXpuGoyiAjFHRK_KeD9HeIW32DcUZrPjNSODXjNo2AA87Bxa8bJXfXMBPBg'}, json={'amount' : 1, 'outcome' : 'NO', 'contractId' : '8b2EC2Au504RG4OCZnqd'}) as resp:
        pokemon = await resp.json(content_type=None)
        return pokemon["betId"]


async def main(length):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for number in range(1, length+1):
            url = 'https://api-nggbo3neva-uc.a.run.app/placebet'
            tasks.append(asyncio.ensure_future(get_pokemon(session, url)))
        original_pokemon = await asyncio.gather(*tasks)
        for pokemon in original_pokemon:
            print(pokemon)

start_time = time.time()
asyncio.run(main(5))
print('--- %s seconds ---' % (time.time() - start_time))

