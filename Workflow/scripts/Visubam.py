from PIL import Image,ImageDraw,ImageFont
import pysam

# create an image
out = Image.new("RGB", (8000, 8000), (255, 255, 255))

# get a font
fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", 40)
# get a drawing context
d = ImageDraw.Draw(out)

# draw multiline text
d.multiline_text((10, 10), "Hello\nWorld", font=fnt, fill=(0, 0, 0))

out.show()






for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=findname(tool)
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)