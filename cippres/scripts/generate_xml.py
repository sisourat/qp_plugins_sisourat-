import sys

na = 1
nb = 1
spin = 1

X = 1
Y = 2
Z = 3

print "<input>"
print ""

print"  <cirun>"
print ""

print "    <general>"
print "      <spin>", spin, "</spin>"
print "      <electron' na='", na, " nb=", nb, " />"
print "    </general>"

print ""

print "    <block>"
print "      <space ne='", X, "' imo='", Y, "' fmo='", Z, "' />"
print "    </block>"

print ""
print "  </cirun>"
print ""

for i in range(5,10):

  print"  <cirun>"

  print "    <general>"
  print "      <spin>", spin, "</spin>"
  print "      <electron' na='", na, " nb=", nb, " />"
  print "    </general>"

  print ""

  print "    <block>"
  print "      <space ne='", 1, "' imo='", i, "' fmo='", i, "' />"
  print "    </block>"

  print "  </cirun>"
  print ""

print "</input>"

