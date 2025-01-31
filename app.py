from dash import Dash, html, dcc, Input, Output
import plotly.express as px
import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np

# Load first dataset
file_path1 = 'data/Relative_PE_UPenn.csv'
df1 = pd.read_csv(file_path1)

# Load second dataset
file_path2 = 'data/Percent_total_UPENN.csv'
df2 = pd.read_csv(file_path2)

# Extract column names (cells) excluding 'target' for both datasets
df1_cells = [col for col in df1.columns if col != 'target']
df2_cells = [col for col in df2.columns if col != 'target']

# Extract unique target values
df1_target_options = df1['target'].unique().tolist()
df2_target_options = df2['target'].unique().tolist()

# Combine unique targets from both datasets
all_targets = list(set(df1_target_options + df2_target_options))

# Use the broader dataset for available cell options
all_cells = list(set(df1_cells + df2_cells))

# Default selections
default_cells = ['mDC', 'pDC']  # Default X-axis (cells)
default_targets = ['CD11b', 'CD20']  # Default Lines (targets)
default_highlight_target = all_targets[0]  # Default highlighted target
default_highlight_cell = all_cells[0]  # Default highlighted cell


# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.MORPH])

# Layout
app.layout = html.Div([
    html.H1('Protein Shake', className='text-center pb-3'),

    # Shared dropdowns for both graphs
    html.Div(className='row', children=[
        html.Div(className='col-6 mx-auto', children=[
            html.Label('Select Targets (Lines for Line Graphs, Dots for Dot Plot):'),
            dcc.Dropdown(
                id='target-selector',
                options=[{'label': target, 'value': target} for target in all_targets], 
                value=default_targets,  
                multi=True,
                className='w-75 mx-auto mb-3'
            ),
        ]),
        html.Div(className='col-6 mx-auto', children=[
            html.Label('Select Cells (X-axis Labels):'),
            dcc.Dropdown(
                id='cell-selector',
                options=[{'label': cell, 'value': cell} for cell in all_cells], 
                value=default_cells,  
                multi=True,
                className='w-75 mx-auto mb-3'
            ),
        ]),
    ]),

    # Row for Line Graphs (Top Row)
    html.Div(className='row', children=[
        html.Div(className='col-6', children=[
            html.H2('Relative PE UPenn', className='text-center'),
            dcc.Graph(id='line-graph-1'),
        ]),
        html.Div(className='col-6', children=[
            html.H2('Percent Total UPenn', className='text-center'),
            dcc.Graph(id='line-graph-2'),
        ]),
    ]),

    # Row for Dot Graphs (Bottom Row)
    html.Div(className='row mt-4', children=[
        # Left bottom dot plot (Select Target)
        html.Div(className='col-6', children=[
            html.Label('Select Target to Highlight in Dot Plot:'),
            dcc.Dropdown(
                id='highlight-target-selector',
                options=[{'label': target, 'value': target} for target in all_targets], 
                value=default_highlight_target,  
                className='w-75 mx-auto mb-3'
            ),
            html.H2('Dot Plot of Targets', className='text-center'),
            dcc.Graph(id='dot-graph-1'),
        ]),

        # Right bottom dot plot (Select Cell)
        html.Div(className='col-6', children=[
            html.Label('Select Cell to Highlight in Dot Plot:'),
            dcc.Dropdown(
                id='highlight-cell-selector',
                options=[{'label': cell, 'value': cell} for cell in all_cells], 
                value=default_highlight_cell,  
                className='w-75 mx-auto mb-3'
            ),
            html.H2('Dot Plot of Cells', className='text-center'),
            dcc.Graph(id='dot-graph-2'),
        ]),
    ]),
])

# Callback to update both line graphs
@app.callback(
    [Output('line-graph-1', 'figure'),
     Output('line-graph-2', 'figure')],
    [Input('target-selector', 'value'),  
     Input('cell-selector', 'value')]   
)
def update_line_graphs(selected_targets, selected_cells):
    if not selected_targets or not selected_cells:
        return px.line(title='Please select at least one target and one cell.'), px.line(title='Please select at least one target and one cell.')

    df1_filtered = df1[df1['target'].isin(selected_targets)][['target'] + selected_cells]
    df1_melted = df1_filtered.melt(id_vars='target', var_name='cell', value_name='expression')
    df1_melted['expression'] = np.arcsinh(df1_melted['expression']/600)

    fig1 = px.line(df1_melted, x='cell', y='expression', color='target', title='Relative PE UPenn')

    df2_filtered = df2[df2['target'].isin(selected_targets)][['target'] + selected_cells]
    df2_melted = df2_filtered.melt(id_vars='target', var_name='cell', value_name='expression')

    fig2 = px.line(df2_melted, x='cell', y='expression', color='target', title='Percent Total UPenn')

    return fig1, fig2

# Callback to update dot plot (Highlight Target)
@app.callback(
    Output('dot-graph-1', 'figure'),
    [Input('highlight-target-selector', 'value')]
)

def update_dot_graph_1(row_name):
    df1_copy = df1.copy(deep=True)
    df1_copy.set_index('target', inplace=True)
    df1_corr = df1_copy.T.corr()
    
    melted_corr1 = df1_corr.reset_index().melt(id_vars='target', var_name='Marker', value_name='CorrelationMFI')
    ans1 = melted_corr1[melted_corr1['target'] == row_name].reset_index(drop=True)

    # compute correlation matrix between all rows: percentages
    df2_copy = df2.copy(deep=True)
    df2_copy.set_index('target', inplace=True)
    df2_corr = df2_copy.T.corr()
    
    melted_corr2 = df2_corr.reset_index().melt(id_vars='target', var_name='Marker', value_name='CorrelationPerc')
    ans2 = melted_corr2[melted_corr2['target'] == row_name].reset_index(drop=True)
    
    ans_both = ans1.merge(ans2, on=['target','Marker'], how='inner')

    fig = px.scatter(
        ans_both,
        x='CorrelationPerc',
        y='CorrelationMFI',
        title=f'Correlation of {row_name} with Markers',
        hover_data=['Marker']
    )
    return fig

# Callback to update dot plot (Highlight Cell)
@app.callback(
    Output('dot-graph-2', 'figure'),
    [Input('highlight-cell-selector', 'value')]
)
def update_dot_graph_2(selected_cell):
    pattern = selected_cell
    
    def custom_summary(group):
        max_mfi_in  = group[group['Celltype'].str.contains(pattern, regex=True)]['MFI'].max()  # Replace with your logic
        max_mfi_out = group[~group['Celltype'].str.contains(pattern, regex=True)]['MFI'].max()  # Replace with your logic
        max_perc_in  = group[group['Celltype'].str.contains(pattern, regex=True)]['Perc'].max()  # Replace with your logic
        max_perc_out  = group[~group['Celltype'].str.contains(pattern, regex=True)]['Perc'].max()  # Replace with your logic

        return pd.Series({
            'Max_MFI_in': max_mfi_in,
            'Max_MFI_out': max_mfi_out,
            'Max_Perc_in': max_perc_in,
            'Max_Perc_out': max_perc_out
        })

    df1_melt = df1.melt(id_vars='target', var_name='Celltype', value_name='MFI')
    df2_melt = df2.melt(id_vars='target', var_name='Celltype', value_name='Perc')
    df_both_melt = df1_melt.merge(df2_melt, on=['target','Celltype'], how='inner')

    df_summary = df_both_melt.groupby('target').apply(custom_summary).reset_index()
    df_filtered = df_summary[(df_summary['Max_MFI_in'] > df_summary['Max_MFI_out']) & (df_summary['Max_Perc_in'] > 0.5)]

    df1_heat = df1[df1['target'].isin(df_filtered['target'])].set_index('target')


    fig = px.imshow(np.arcsinh(df1_heat/600), 
                    labels=dict(x="Features", y="Categories", color="Value"),
                    color_continuous_scale='viridis',  # Change color scale
                    text_auto=True)  # Show values inside cells

    return fig


# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
