import React, { Component } from 'react'

export class GeneralStudy extends Component {
    constructor(props) {
        super(props)

        this.state = {

        }
    }


    render() {
        return (
            <div className="container-fluid">
                <div className="row">
                    <div className="col-lg-12 col-md-12 col-sm-12">
                        
                        <div class="card">
                        <b className="investigation">Investigation</b>
                            <div class="card-body">
                                <div id="main">
                                    <div class="container">
                                        <div class="accordion" id="faq">
                                            <div class="card">
                                                <div class="card-header" id="faqhead1">
                                                    <a href="#" class="btn btn-header-link" data-toggle="collapse" data-target="#faq1"
                                                        aria-expanded="true" aria-controls="faq1"><b>General Information</b></a>
                                                </div>

                                                <div id="faq1" class="collapse show" aria-labelledby="faqhead1" data-parent="#faq">
                                                    <div class="card-body">
                                                        <form>
                                                            <div className="form-group">
                                                                <label for="email">Title</label>
                                                                <input type="text" class="form-control" placeholder="Enter Title" />
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="comment">Description</label>
                                                                <textarea class="form-control" placeholder="Write Description" rows="3"></textarea>
                                                            </div>
                                                            <div className="d-flex ">
                                                                <div class="form-group">
                                                                    <label for="Submission Date">Submission Date</label>
                                                                    <input type="date" class="form-control" id="Submission Date" />
                                                                </div>
                                                                <div class="form-group ml-4">
                                                                    <label for="Release Date">Public Release Date</label>
                                                                    <input type="date" class="form-control" id="Release  Date" />
                                                                </div>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Publications">Publications</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Publications">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Document Licenses">ISA Document Licenses</label>
                                                                <a type="button" className="btn float-right" href="#"><i class="fa fa-plus"></i>AddNew</a>
                                                                <select class="form-control" id="Document Licenses">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Contacts">Contacts</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Contacts">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>

                                                        </form>
                                                    </div>
                                                </div>
                                            </div>
                                            <div class="card">
                                                <div class="card-header" id="faqhead2">
                                                    <a href="#" class="btn btn-header-link collapsed" data-toggle="collapse" data-target="#faq2"
                                                        aria-expanded="true" aria-controls="faq2"><b>Study 1</b></a>
                                                </div>

                                                <div id="faq2" class="collapse" aria-labelledby="faqhead2" data-parent="#faq">
                                                    <div class="card-body">
                                                        <form>
                                                            <div className="form-group">
                                                                <label for="email">Title</label>
                                                                <input type="text" class="form-control" placeholder="Enter Title" />
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="comment">Description</label>
                                                                <textarea class="form-control" placeholder="Write Description" rows="3"></textarea>
                                                            </div>

                                                            <div class="form-group">
                                                                <label for="Publications">Publications</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Publications">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            {/* <div class="form-group">
                                                                <label for="Document Licenses">ISA Document Licenses</label>
                                                                <a type="button" className="btn float-right" href="#"><i class="fa fa-plus"></i>AddNew</a>
                                                                <select class="form-control" id="Document Licenses">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div> */}
                                                            <div class="form-group">
                                                                <label for="Contacts">Contacts</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Contacts">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <label for="Study Type">Study Type</label>
                                                            <div className="d-flex ">

                                                                <button type="button" className="btn  btn-light">Observational Study</button>
                                                                <button type="button" className="btn  btn-light ml-3">Intervantional Study</button>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Publications">Design Type</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Publications" placeholder="Select Design Type">
                                                                    <option>Select Design Type</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Publications">Factor Name</label>

                                                                <input type="text" class="form-control" placeholder="Enter Factor Name" />

                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Publications">Factor Type</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Publications">
                                                                    <option>Select Factor Type</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <div class="form-group">
                                                                <label for="Publications">Protocols</label>
                                                                <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                <select class="form-control" id="Publications">
                                                                    <option>1</option>
                                                                    <option>2</option>
                                                                    <option>3</option>
                                                                    <option>4</option>
                                                                </select>
                                                            </div>
                                                            <div class="accordion" id="faqTwo">
                                                                <div class="card">
                                                                    <div class="card-header" id="faqhead3">
                                                                        <a href="#" class="btn btn-header-link collapsed" data-toggle="collapse" data-target="#faqTwo"
                                                                            aria-expanded="true" aria-controls="faq3"><b>Assay 1</b></a>
                                                                    </div>

                                                                    <div id="faqTwo" class="collapse" aria-labelledby="faqhead3" data-parent="#faqTwo">
                                                                        <div class="card-body">
                                                                            <form>
                                                                                <div className="form-group">
                                                                                    <label for="email">Title</label>
                                                                                    <input type="text" class="form-control" placeholder="Enter Title" />
                                                                                </div>

                                                                                <div className="d-flex ">
                                                                                    <div class="form-group">
                                                                                        <label for="Start Date">Start Date</label>
                                                                                        <input type="date" class="form-control" id="Start Date" />
                                                                                    </div>

                                                                                    <div class="form-group ml-4">
                                                                                        <label for="End Date">End Date</label>
                                                                                        <input type="date" class="form-control" id="End Date" />
                                                                                    </div>
                                                                                    <div class="form-group ml-4">
                                                                                        <label for="Run Order">Run Order</label>
                                                                                        <input type="text" class="form-control" id="Run Order" />
                                                                                    </div>
                                                                                </div>
                                                                                <div class="form-group">
                                                                                    <label for="Performer">Performer</label>
                                                                                    <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                                    <select class="form-control" id="Performer">
                                                                                        <option>Select Performer</option>
                                                                                        <option>2</option>
                                                                                        <option>3</option>
                                                                                        <option>4</option>
                                                                                    </select>
                                                                                </div>
                                                                                <div class="card mb-4">
                                                                                    <div class="card-body"><label>Inputs</label>
                                                                                        <div className="row mb-4">
                                                                                            <div className="col-lg-3">
                                                                                                <label for="Bio">Digital Object</label>
                                                                                                <button type="button" className="btn btn-light"><i class="fa fa-check-circle"></i>File Uploaded</button>
                                                                                            </div>
                                                                                            <div className="col-lg-6">
                                                                                                <div className="form-group">
                                                                                                    <label for="Bio">Bio Material Object</label>
                                                                                                    <input type="text" class="form-control" placeholder="Bio Material Object" />
                                                                                                </div>
                                                                                            </div>
                                                                                        </div>

                                                                                        <div className="d-flex justify-content-center pb-4">
                                                                                            <button type="button" className="btn inputButton"><i class="fa fa-plus-circle fa_modified"></i>Add New Inputs</button>
                                                                                        </div>
                                                                                    </div>

                                                                                </div>
                                                                                <div class="card">
                                                                                    <div class="card-body"><label>Outputs</label>
                                                                                        <div className="row mb-4">
                                                                                            <div className="col-lg-3">
                                                                                                <label for="Bio">Digital Object</label>
                                                                                                <button type="button" className="btn btn-light"><i class="fa fa-file-upload"></i>File Upload</button>
                                                                                            </div>
                                                                                            <div className="col-lg-6">
                                                                                                <div className="form-group">
                                                                                                    <label for="Bio">Bio Material Object</label>
                                                                                                    <input type="text" class="form-control" placeholder="Bio Material Object" />
                                                                                                </div>
                                                                                            </div>
                                                                                        </div>

                                                                                        <div className="d-flex justify-content-center pb-4">
                                                                                            <button type="button" className="btn inputButton"><i class="fa fa-plus-circle fa_modified"></i>Add New Outputs</button>
                                                                                        </div>
                                                                                    </div>

                                                                                </div>
                                                                                <label className="mt-3">Measurement Type</label>
                                                                                <div className="form-group">
                                                                                    <label for="Bio">Name</label>
                                                                                    <input type="text" class="form-control" placeholder="Enter Name" />
                                                                                </div>
                                                                                <div class="form-group">
                                                                                    <label for="Performer">Annotation</label>
                                                                                    <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                                    <select class="form-control" id="Performer">
                                                                                        <option>Select Annotation</option>
                                                                                        <option>2</option>
                                                                                        <option>3</option>
                                                                                        <option>4</option>
                                                                                    </select>
                                                                                </div>

                                                                                <label>Technology type</label>
                                                                                <div className="form-group">
                                                                                    <label for="Bio">Name</label>
                                                                                    <input type="text" class="form-control" placeholder="Enter Name" />
                                                                                </div>
                                                                                <div class="form-group">
                                                                                    <label for="Performer">Annotation</label>
                                                                                    <button type="button" className="btn float-right"><i class="fa fa-plus"></i>AddNew</button>
                                                                                    <select class="form-control mb-5" id="Performer">
                                                                                        <option>Select Annotation</option>
                                                                                        <option>2</option>
                                                                                        <option>3</option>
                                                                                        <option>4</option>
                                                                                    </select>
                                                                                </div>
                                                                                <div className="d-flex justify-content-center pb-4">
                                                                                            <button type="button" className="btn inputButton"><i class="fa fa-plus-circle fa_modified"></i>Add New Assay</button>
                                                                                        </div>
                                                                            </form>


                                                                        </div>
                                                                    </div>
                                                                </div>
                                                            </div>
                                                        </form>

                                                    </div>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        )
    }
}

export default GeneralStudy


