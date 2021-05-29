
import React, { Component } from 'react'

export class LeftSidebar extends Component {
    constructor(props) {
        super(props)
    
        this.state = {
             
        }
    }
    
    render() {
        return (
<div className="sideBar"> 
      
<nav class="navbar">
<ul className="navbar-nav pt-3">
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-file-code leftSidebarIcon"></i>Dashboard</a>
  </li>
  <li className="nav-item">
    <a class="nav-link" href="#"><i class="fa fa-file-code leftSidebarIcon"></i>Source code</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-lightbulb leftSidebarIcon"></i>Ideas</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-users leftSidebarIcon"></i>Contacts</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-user leftSidebarIcon"></i>Agents</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-newspaper leftSidebarIcon"></i>Articles</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-cog leftSidebarIcon"></i>Settings</a>
  </li>
  <li className="nav-item">
    <a className="nav-link" href="#"><i class="fa fa-ring leftSidebarIcon"></i>Subscriptions</a>
  </li>
</ul>

</nav>
<div className="sidebarBox">
<span><i class="fa fa-file-alt sidebarboxIcon"></i></span>
  <div className="sideboxContent">
  
  <span>Need Help?</span>
  <p>Please Check our docs</p>
<button type="button" className="btn leftBtn">Documentation</button>                                          
</div>
</div>
</div>    
        )
    }
}

export default LeftSidebar
