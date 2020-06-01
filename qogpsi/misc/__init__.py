from random import choice
from string import ascii_lowercase
# Misc functions

def generate_string(n=10):
	return "".join(choice(ascii_lowercase) for i in range(n))
