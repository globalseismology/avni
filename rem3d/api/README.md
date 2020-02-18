# Using the API

## initial setup
1. get an api key. Log into https://maurya.umd.edu/ and visit https://maurya.umd.edu/tools/applets/ to get an api key
2. (optional) If not on the main branch containing the api, you can check out the api branch from the rem3d repo. Depending on how you've installed the rem3d package, you may need to reinstall after switching branches (e.g., `python setup.py install --user` from the conda environment you're using for rem3d).
3. store your api key:
 * open up a python terminal
 * import the api client:
 ```
 from rem3d.api.client import Client as r3d
 ```
 * initialize a client connection object with
  ```
  conn=r3d(api_key='API_KEY')
  ```
  where ``'API_KEY'`` is your api keystring
 * store the api key in the default api config file by calling
 ```
 conn.setApiConfig()
 ```

At this point, your api key will be stored in the config file, `rem3d/config/api.ini` and future calls to initialize a client will pull the key from this file (i.e., you can do `conn=r3d()` to initialize).

After completing step 3 above, you can run `conn.checkUserStats()` to validate that your key works.

## using the api

The basic procedure for using the API is to:
1. initialize the client connection object
2. initialize api classes with the connection object
3. call methods within api classes using dictionary arguments. All the methods will return results in a dictionary, with descriptive keys.

Some notes:
- Each call the api is tracked, and you can call `conn.checkUserStats()` at any point to see your call count. The current max calls per 30 days is set at 500. If you exceed this count you'll get a warning, and subsequent calls will lock your account for 15 days. At present we do not limit the rate of calls, but please be considerate of the shared resource and consider self-limiting calls. 
- Examples scripts using the api are stored in `examples/Scripts/api_*`
