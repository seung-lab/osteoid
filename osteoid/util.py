from .skeleton import Skeleton

def load(filename:str) -> Skeleton:
	with open(filename, "rb") as f:
		data = f.read()

	if filename.endswith("swc"):
		return Skeleton.from_swc(data.decode("utf8"))
	elif filename.endswith("ostd"):
		return Skeleton.from_ostd(data)
	else:
		return Skeleton.from_precomputed(data)

def save(filename:str, skel:Skeleton):
	if filename.endswith("swc"):
		binary = skel.to_swc()
	elif filename.endswith("ostd"):
		binary = skel.to_ostd()
	else:
		binary = skel.to_precomputed()

	with open(filename, "wb") as f:
		f.write(binary)
