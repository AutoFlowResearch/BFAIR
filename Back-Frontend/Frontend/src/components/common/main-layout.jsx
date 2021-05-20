
import React, { Component } from 'react'
import SidebarArea from './sidebarArea/index';
import ContentArea from './contentArea/index';
import GeneralinvestigationSidebar from '../pages/investigation-form-layout/investigation-form-sidebar'
import GeneralStudy from '../pages/investigation-form-layout/general-form-layout'

export class MainLayout extends Component {
    constructor(props) {
        super(props)
    
        this.state = {
             
        }
    }
    
    render() {
        return (
            
            <div className="container-fluid p-0">
                <div className="card bg-container Maincard">
                    <div className="row m-0 mt-4" style = {{paddingBottom: '18px'}}>
                    <div className=" col-lg-2 col-md-2 col-sm-2">
                            <div className="brand">
                                <div className="logo"><b>B</b></div>
                                <span><b>BFAIR</b></span>
                            </div>
                        </div>
                        <div className="col-lg-10 col-md-10 col-sm-10 ">
                        <div class="d-flex justify-content-between">
                            <h5 className="font-weight-bold">Dashboard</h5>
                            {/* <div>
                            <p>lawrence</p>
                            </div> */}
                            <div className="brand">
                                <div className="headerUSer"><b>B</b></div>
                                <span>lawrence</span>
                            </div>
                        </div>
                        </div>
                    </div>
                <div className="row m-0 ">
                <div className=" p-0 col-lg-2 col-md-2 col-sm-2">
                <SidebarArea></SidebarArea>
                
                </div>
                <div className="col-lg-10 col-md-10 col-sm-10 p-0">
                 
                    <ContentArea></ContentArea>
                </div>
            </div>
            {/* <div className="row m-0">
                    <div className=" p-0 col-lg-2 col-md-2 col-sm-2">
                    <GeneralinvestigationSidebar></GeneralinvestigationSidebar>
                    </div>
                    <div className="col-lg-10 col-md-10 col-sm-10 p-0">
                        <GeneralStudy></GeneralStudy>
                    </div>
                </div> */}
                </div>
            </div>
        )
    }
}

export default MainLayout
