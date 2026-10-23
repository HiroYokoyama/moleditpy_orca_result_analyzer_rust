# Shared maintenance tools

Do not edit copied scripts independently. Their tests run in the canonical repositories, not in every consumer CI matrix.

API checker: the canonical engine remains in [moleditpy-plugins/api-checker](https://github.com/HiroYokoyama/moleditpy-plugins/tree/main/api-checker), version 1.1.3. `api-checker.lock.json` pins its source commit and hash. Keep this plugin's API compatibility check; engine unit tests run only in the collection.
