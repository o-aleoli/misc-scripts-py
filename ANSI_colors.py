#!/usr/bin/python3
colors = list(map(str, range(16, 232)))
text = "color!"
ctrl_seq_intro = "\u001b["
ctrl_seq_reset_all = "\u001b[0m"
fg_code = "38;5;"
bg_code = "48;5;"

for i in range(len(colors)):
    print(
        ctrl_seq_intro\
        + fg_code + colors[i] + ";"\
        + bg_code + colors[-i]\
        + "m" + text\
        + ctrl_seq_reset_all
    )
