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

pattern = 'B cell'  # Change this pattern as needed
df_summary = df_both_melt.groupby('target').apply(custom_summary).reset_index()
df_filtered = df_summary[(df_summary['Max_MFI_in'] > df_summary['Max_MFI_out']) & (df_summary['Max_Perc_in'] > 0.5)]

df1_heat = df1[df1['target'].isin(df_filtered['target'])].melt(id_vars='target', var_name='Celltype', value_name='MFI')
df1_heat['MFI'] = np.arcsinh(df1_heat['MFI']/600)
print(df1_heat)

fig = px.imshow(df1_heat, labels=dict(x="Celltype", y="target", color="MFI"),
                color_continuous_scale='viridis', text_auto=True)

# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.MORPH])

# Layout
app.layout = html.Div([
    html.H1('Comparison of Two Datasets', className='text-center pb-3'),

    # Shared dropdowns for both graphs
    html.Div(className='row', children=[
        html.Div(className='col-6 mx-auto', children=[
            html.Label('Select Targets (Lines for Line Graphs, Dots for Dot Plot):'),
            dcc.Dropdown(
                id='target-selector',
                options=[{'label': target, 'value': target} for target in all_targets], 
                value=default_targets,  # Default selected targets
                multi=True,
                className='w-75 mx-auto mb-3'
            ),
            html.Label('Select Cells (X-axis Labels):'),
            dcc.Dropdown(
                id='cell-selector',
                options=[{'label': cell, 'value': cell} for cell in all_cells], 
                value=default_cells,  # Default selected values
                multi=True,
                className='w-75 mx-auto mb-3'
            ),
        ]),
    ]),

    html.Div(className='row', children=[
        # Left graph section (Relative_PE_UPenn)
        html.Div(className='col-6', children=[
            html.H2('Relative PE UPenn', className='text-center'),
            dcc.Graph(id='line-graph-1'),
        ]),

        # Right graph section (Percent_total_UPENN)
        html.Div(className='col-6', children=[
            html.H2('Percent Total UPenn', className='text-center'),
            dcc.Graph(id='line-graph-2'),
        ]),
    ]),

    # Dot plot section
    html.Div(className='row mt-4', children=[
        html.Div(className='col-6 mx-auto', children=[
            html.Label('Select Target to Highlight in Dot Plot:'),
            dcc.Dropdown(
                id='highlight-target-selector',
                options=[{'label': target, 'value': target} for target in all_targets], 
                value=default_highlight_target,  # Default highlighted target
                className='w-75 mx-auto mb-3'
            ),
        ]),
    ]),

    html.Div(className='row', children=[
        html.Div(className='col-12', children=[
            html.H2('Dot Plot of Targets', className='text-center'),
            dcc.Graph(id='dot-graph'),
        ]),
    ]),
])

# Callback to update both line graphs
@app.callback(
    [Output('line-graph-1', 'figure'),
     Output('line-graph-2', 'figure')],
    [Input('target-selector', 'value'),  # Selected targets (Lines)
     Input('cell-selector', 'value')]   # Selected cells (X-axis)
)
def update_line_graphs(selected_targets, selected_cells):
    """
    Updates both line graphs based on shared target and cell selections.

    Args:
        selected_targets (list): List of selected line values.
        selected_cells (list): List of selected x-axis values.

    Returns:
        tuple: Updated line charts for both datasets.
    """
    if not selected_targets or not selected_cells:
        return px.line(title='Please select at least one target and one cell.'), px.line(title='Please select at least one target and one cell.')

    # Process the first dataset
    df1_filtered = df1[df1['target'].isin(selected_targets)][['target'] + selected_cells]
    df1_melted = df1_filtered.melt(id_vars='target', var_name='cell', value_name='expression')
    df1_melted['expression'] = np.arcsinh(df1_melted['expression']/600)

    fig1 = px.line(df1_melted, x='cell', y='expression', color='target', title='Relative PE UPenn')
    fig1.update_xaxes(title_text='Selected Cells (X-axis)')
    fig1.update_yaxes(title_text='Expression Levels')

    # Process the second dataset
    df2_filtered = df2[df2['target'].isin(selected_targets)][['target'] + selected_cells]
    df2_melted = df2_filtered.melt(id_vars='target', var_name='cell', value_name='expression')

    fig2 = px.line(df2_melted, x='cell', y='expression', color='target', title='Percent Total UPenn')
    fig2.update_xaxes(title_text='Selected Cells (X-axis)')
    fig2.update_yaxes(title_text='Expression Levels')

    return fig1, fig2

# Callback to update dot plot with highlighted target
@app.callback(
    Output('dot-graph', 'figure'),
    [Input('highlight-target-selector', 'value')]  # Selected target to highlight
)
def update_dot_graph(row_name):
    """
    Updates the dot graph to highlight a selected target in pink.

    Args:
        highlight_target (str): Selected target to highlight.

    Returns:
        plotly.graph_objs._figure.Figure: Updated dot graph.
    """
    # Combine both datasets
    # df_combined = pd.concat([df1, df2])

    # Add a column to distinguish the highlighted target
    # df_combined['highlight'] = df_combined['target'].apply(lambda x: 'Highlighted' if x == highlight_target else 'Other')

    # compute correlation matrix between all rows: MFI
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

    # Create dot plot
    fig = px.scatter(
        ans_both,
        x='CorrelationPerc',
        y='CorrelationMFI',
        title='Correlation between markers and ' + row_name,
        hover_data=['Marker']
    )

    # fig.update_xaxes(title_text='Targets')
    # fig.update_yaxes(title_text='Expression Levels')

    return fig


# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
