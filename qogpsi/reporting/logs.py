
def log_title(title="Insert A Title Here", **kwargs):
	indent = kwargs.get('indent', "")
	if kwargs.get('new_line', False):
		print("\n")
	print("%s-------- %s --------" % (indent, title))

def log_subtitle(subtitle="Insert A Subtitle Here", **kwargs):
	indent = kwargs.get('indent', "")
	if kwargs.get('new_line', False):
		print("\n")
	print("%s    ---- %s ----" % (indent, subtitle))

def log_dictionary(dictionary={}, **kwargs):
	indent = kwargs.get('indent', "")
	if kwargs.get('new_line', False):
		print("\n")
	if (dictionary != {}):
		for k,v in dictionary.items():
			print("%s\t%s : %s" % (indent, k, v))

def log_text(text, **kwargs):
	indent = kwargs.get('indent', "")
	if kwargs.get('new_line', False):
		print("\n")
	print("%s%s" % (indent, text))