import React from 'react';
import { Link } from 'react-router-dom';
let baseUrl = 'http://127.0.0.1:5000';

class SampleSubmission extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            sampleModule: '',
            sampleSelect: '',
            sampleData: [],
            investigationData : [],
            typeName: '',
            method: '',
            description: '',
            research_center: '',
            investigationTitleInput: '',
            investigationDescInput: '',
            investigationDescSelect: '',
            investigation_Title: '',
            investigation_Id: ''
        };

        this.sampleChange = this.sampleChange.bind(this);
        this.investigationTitleChange = this.investigationTitleChange.bind(this)
        this.deleteSample = this.deleteSample.bind(this);
        this.deleteInvestigation = this.deleteInvestigation.bind(this);
        this.clearSample = this.clearSample.bind(this);
        this.changeText = this.changeText.bind(this);
        this.typeChange = this.typeChange.bind(this);
        this.methodChange = this.methodChange.bind(this);
        this.centerChange = this.centerChange.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleSampleModuleSubmit = this.handleSampleModuleSubmit.bind(this);
        this.clearInvestigatingForm = this.clearInvestigatingForm.bind(this);
        this.handleInvestigationSubmit = this.handleInvestigationSubmit.bind(this);
        this.updateInvestigation = this.updateInvestigation.bind(this);
        
    }

    componentDidMount() {
        this.getSampleDataList();
        this.getInvestigationList();
    }

    clearSample() {
        this.setState({
            sampleModule: ''
        });
        console.log(this.state.sampleModule);
    }

    clearInvestigatingForm() {
        this.setState({
            investigationTitleInput: ''
        });
    }

    sampleChange(event) {
        this.setState({
            sampleSelect: event.target.value
        });
        console.log(event.target.value);
    }

    investigationTitleChange(event){
        this.setState({
            investigationSelect: event.target.value
        });
        let selectedObj = this.state.investigationData.filter(function (el) {
            return el.Investigation_id == event.target.value
          });
        this.setState({
            investigationDescSelect: selectedObj[0].Description,
            investigation_Id: selectedObj[0].Investigation_id,
            investigation_Title: selectedObj[0].Title
        });

    }

    deleteSample(event) {
        console.log('sample selete', this.state.sampleSelect);
        const requestOptions = {
            method: 'DELETE',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ Sample_id: this.state.sampleSelect })
        };
        fetch(`${baseUrl}/Sample_api`, requestOptions)
            .then((response) => {
                alert('Sample successfully deleted')
                this.getSampleDataList();
                response.json()
            } )
    }

    deleteInvestigation(event) {
        // console.log('sample selete', this.state);
        const requestOptions = {
            method: 'DELETE',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ Investigation_id: this.state.investigationSelect })
        };
        fetch(`${baseUrl}/investigation_api`, requestOptions)
            .then((response) => {
                alert('Investigation record successfully deleted')
                this.getInvestigationList();
                response.json()
            } )
    }


    typeChange(event) {
        this.setState({
            typeName: event.target.value
        });
        console.log(event.target.value);
    }

    methodChange(event) {
        this.setState({
            method: event.target.value,
            description: 'auto_populated'
        });
    }

    centerChange(event) {
        this.setState({
            method: event.target.value
        });
        console.log(event.target.value);
    }

    handleSubmit(e) {
        e.preventDefault();
        const { data } = this.state;
        let obj = {
            typeName: this.state.typeName,
            method: this.state.method,
            description: this.state.description,
            research_center: this.state.research_center
        }
    }

    updateInvestigation(){
        console.log(this.state);
        // e.preventDefault();
        const requestOptions = {
            method: 'PUT',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ Title: this.state.investigation_Title,Description: this.state.investigationDescInput, Investigation_id: this.state.investigation_Id })
        };
        fetch(`${baseUrl}/investigation_api`, requestOptions)
            .then(response => response.json())
            .then((data) => {
                alert('Investigation record updated successfully')
                this.getInvestigationList();
            })

    }

    handleInvestigationSubmit(e){
        e.preventDefault();
        const requestOptions = {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ Title: this.state.investigationTitleInput,Description: this.state.investigationDescInput })
        };
        fetch(`${baseUrl}/investigation_api`, requestOptions)
            .then(response => response.json())
            .then((data) => {
                alert(data.message)
                this.getInvestigationList();
            })
    }

    setSampleModuleInput(e) {
        this.setState({
            sampleModule: e
        });
    }

    investigationSetTitleInput(e){
        this.setState({
            investigationTitleInput: e
        });
    }

    investigationSetDescInput(e){
        this.setState({
            investigationDescInput: e
        });
    }

    handleSampleModuleSubmit(e) {
        e.preventDefault();
        console.log('sample name', this.state.sampleModule);
        const requestOptions = {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ Sample_name: this.state.sampleModule })
        };
        fetch(`${baseUrl}/Sample_api`, requestOptions)
            .then(response => response.json())
            .then((data) => {
                alert(data.message)
                this.getSampleDataList();
            })
    }

    getSampleDataList(){
        fetch(`${baseUrl}/Sample_api`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
             this.setState({
                 sampleData: data.Sample
             })
        })
    }

    getInvestigationList(){
        fetch(`${baseUrl}/investigation_api`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
             this.setState({
                investigationData: data.Investigation
             })
        })

    }

    changeText(event){
        this.setState(
            {investigationDescInput : event.target.value}
        );
    }

    render() {
        let optionTemplate = this.state.sampleData.map(v => (
            <option value={v.id}>{v.Sample_name}</option>
          ));

        let optionInvestigation = this.state.investigationData.map(v => (
            <option value={v.Investigation_id}>{v.Title}</option>
          ));
        return (
            <div className="container">
                <ul className="nav nav-tabs" role="tablist">
                    <li className="nav-item">
                        <a className="nav-link active" data-toggle="tab" href="#SampleModule">Sample Module</a>
                    </li>
                    <li className="nav-item">
                        <a className="nav-link" data-toggle="tab" href="#investigation">Investigation</a>
                    </li>
                    <li className="nav-item">
                        <a className="nav-link" data-toggle="tab" href="#study">Study Module</a>
                    </li>
                    <li className="nav-item">
                        <a className="nav-link" data-toggle="tab" href="#protocol">Protocols Module</a>
                    </li>
                    <li className="nav-item">
                        <a className="nav-link" data-toggle="tab" href="#assay">Assay Module</a>
                    </li>


                </ul>

                {/* Sample Module Form */}
                <div className="tab-content">
                    <div id="SampleModule" className="container tab-pane active">
                        <div className="col-md-4 col-md-offset-3">
                            <h2>Sample Module</h2>
                            <form name="form" onSubmit={this.handleSampleModuleSubmit}>
                                <div className={'form-group'}>
                                    <label>Sample</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>

                                <div className="form-group">
                                    <button className="btn btn-primary">Add</button>
                                    <button type="button" onClick={this.clearSample} className="btn btn-light">Cancel</button>
                                </div>
                            </form>

                            <div className={'form-group'} >
                                <label htmlFor="sample" style={{ marginTop: '40px' }}>Sample List</label>
                                <select className="form-control" value={this.state.value}  onChange={this.sampleChange}>
                                <option disabled selected>Select Sample</option>
                                {optionTemplate}
                                </select>
                            </div>
                            <div className="form-group">
                                <button onClick={this.deleteSample} className="btn btn-danger">Delete</button>
                            </div>
                        </div>
                    </div>

                    {/* Sample Submission */}

                    <div id="investigation" className="container tab-pane fade">

                        <div className="col-md-4 col-md-offset-3">
                            <h2>Investigation</h2>
                            <form name="form" onSubmit={this.handleInvestigationSubmit}>
                            <div className={'form-group'}>
                                    <label>Title</label>
                                    <input type="text" onInput={e => this.investigationSetTitleInput(e.target.value)} className="form-control" name="investigationTitleInput" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Description</label>
                                    <input type="text" onInput={e => this.investigationSetDescInput(e.target.value)} className="form-control" name="investigationDescInput" />
                                </div>
                                <div className="form-group">
                                    <button className="btn btn-primary">Add</button>
                                    <button type="button" onClick={this.clearInvestigatingForm} className="btn btn-light">Cancel</button>
                                </div>
                                <div className={'form-group'} >
                                <label htmlFor="sample" style={{ marginTop: '40px' }}>Title List</label>
                                <select className="form-control" value={this.state.value}  onChange={this.investigationTitleChange}>
                                <option disabled selected>Select Title</option>
                                {optionInvestigation}
                                </select>
                            </div>
                                <div className={'form-group'}>
                                    <label htmlFor="password">Description</label>
                                    <input type="text" value={this.state.investigationDescSelect} onChange = {this.changeText} className="form-control" name="description" />
                                </div>
                                <div className="form-group">
                                    <button className="btn btn-primary" type="button" onClick={this.updateInvestigation}>Update</button>
                                    <button type="button" onClick={this.deleteInvestigation} className="btn btn-danger">Delete</button>
                                </div>
                            </form>
                        </div>
                    </div>
                    <div id="study" className="container tab-pane fade">
                    <div className="col-md-4 col-md-offset-3">
                            <h2>Study</h2>
                            <form name="form" onSubmit={this.handleSampleModuleSubmit}>
                                <div className={'form-group'}>
                                    <label>Title</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Description</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Submission Date</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Release Date</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Full Name</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Assay</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Contact</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Design type</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Factors</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Protocols</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                            </form>
                            <div className="form-group">
                                <button className="btn btn-primary">Save</button>
                            </div>
                        </div>
                        
                    </div>
                    <div id="protocol" className="container tab-pane fade">
                     <div className="col-md-4 col-md-offset-3">
                            <h2>Protocol</h2>
                            <form name="form" onSubmit={this.handleSampleModuleSubmit}>
                                <div className={'form-group'}>
                                    <label>Title</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Type</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Description</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>URI</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Version</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Parameter</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Component Name</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Component Type</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                              
                            </form>
                            <div className="form-group">
                                <button className="btn btn-primary">Save</button>
                            </div>
                        </div>
                        
                    </div>
                    <div id="assay" className="container tab-pane fade">
                     <div className="col-md-4 col-md-offset-3">
                            <h2>Assay</h2>
                            <form name="form" onSubmit={this.handleSampleModuleSubmit}>
                                <div className={'form-group'}>
                                    <label>Title</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Measurement Type</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Technology Type</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                                <div className={'form-group'}>
                                    <label>Technology Platform</label>
                                    <input type="text" onInput={e => this.setSampleModuleInput(e.target.value)} className="form-control" name="description" />
                                </div>
                              
                            </form>
                            <div className="form-group">
                                <button className="btn btn-primary">Save</button>
                            </div>
                        </div>
                        
                    </div>
                </div>
            </div>
        );
    }
}

export default SampleSubmission;