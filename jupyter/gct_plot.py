import plotly.express as px
import plotly.offline as po
import numpy as np
from scgid.scripts.infotable import InfoTable
it = InfoTable()
it.load("/home/aimzez/development/scgid/test_data/stylopage_41_scgid_output/gct/stylopage_41.infotable.tsv")
it.df = it.df.assign(scaled_evalue = np.log(it.df.evalue))

fig = px.scatter(x=np.log(it.df.coverage),
                 y=it.df.gc,
                 color=it.df.pertinent_taxlvl,
                 size=abs(it.df.scaled_evalue),
                 opacity=0.7
                )
fig.update_layout(
    yaxis_title="GC Content",
    xaxis_title="log(Coverage)"
)
po.plot(fig, filename="your_gct_plot.html", auto_open = True)

