<h1>Mint Simulation</h1>
Python code to create hypothetical die charts and analyse their properties.

<h2>Usage</h2>
The code can be run from the command line. The required arguments are the number of workstations (W), the number of dies in the obverse die box (odb), and the number of dies in the reverse die box (rdb) (in that order). A file name for the output must also be specified. E.g.

```
python3 mint_simulation.py 1 1 1 --output file_name.xlsx
```
The following arguments are optional:
<table>
  <tr>
    <td>-k</td>
    <td>The shape of the gamma distribution (Default: k=1 suggested by Esty 2011. k=2 suggested by Carter 1983)</td>
  </tr>
  <tr>
    <td>-theta</td>
    <td>Average lifetime of reverse dies in minutes (Default: theta=450)</td>
  </tr>
  <tr>
    <td>-a</td>
    <td>Scale factor, how much longer obverses last than reverses (Default: a=1.5)</td>
  </tr>
  <tr>
    <td>-t</td>
    <td>Length of work period in minutes (Default: t=600)</td>
  </tr>
  <tr>
    <td>-D</td>
    <td>Maximum number of obverse dies to be created (Default: D=20)</td>
  </tr>
  <tr>
    <td>-o_loss</td>
    <td>Percentage of obverse dies to delete from full chart (Default: o_loss=0)</td>
  </tr>
  <tr>
    <td>-r_loss</td>
    <td>Percentage of reverse dies to delete from full chart (Default: r_loss=0)</td>
  </tr>
  <tr>
    <td>-e_loss</td>
    <td>Percentage of edges to delete from full chart (Default: e_loss=0)</td>
  </tr>
  <tr>
    <td>-i</td>
    <td>Number of iterations (Default: i=10000)</td>
  </tr>
</table>

Example command line with optional arguments:
<br>
```
python3 mint_simulation.py 1 2 2 -a 3 -t 400 -D 10 --output file_name.xlsx
```
