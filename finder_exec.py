# file to execute code from terminal

import finder_code as code

def help():
    print("")
    print("Here is the list of commands:")
    print("  help : shows this page")
    print("  quit : quits the program")
    print("")

def search():
    pass


while True:
    command = input("c: ")
    if command == "q" or command == "quit":
        break
    elif command == "help":
        help()
    elif command == "search":
        search()
    elif command == "f1":
        code.func()
    else:
        print("Unknown command... type help to see list.")
    
