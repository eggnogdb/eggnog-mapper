import time
from bottle import Bottle, run, get, post, request, route, response, abort, hook, error, HTTPResponse

trees = {}

@route('test')
def test_change():
    trees[time.ctime()] = time.time()
    print trees

if __name__ == "__main__":
    run(host="localhost", port=8989, server='cherrypy')
