from random import choice
from string import ascii_lowercase

def generate_string(n=10):
	return "".join(choice(ascii_lowercase) for i in range(n))

class Option(object):
	def __init__(self, name):
		self.name = name
		self.id = generate_string()